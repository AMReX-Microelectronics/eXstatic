#include "Transport.H"

//#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
//#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/SelectWarpXUtils/WarpXConst.H"
#include "../../Utils/CodeUtils/CodeUtil.H"
#include "../../Code.H"
#include "../../Input/GeometryProperties/GeometryProperties.H"
//#include "Input/MacroscopicProperties/MacroscopicProperties.H"
//#include "Input/BoundaryConditions/BoundaryConditions.H"
//
#include <AMReX.H>
//#include <AMReX_ParmParse.H>
//#include <AMReX_Parser.H>
//#include <AMReX_RealBox.H>
//#include <AMReX_MultiFab.H>
////#include <AMReX_MultiFabUtil.H>


using namespace amrex;

c_TransportSolver::c_TransportSolver()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_TransportSolver() Constructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    ReadData();
  
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_TransportSolver() Constructor************************\n";
#endif
}


c_TransportSolver::~c_TransportSolver()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_TransportSolver() Destructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    vp_CNT.clear();
    vp_Graphene.clear();
    //vp_Silicon.clear();
  
#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_TransportSolver() Destructor************************\n";
#endif
}

void 
c_TransportSolver::ReadData() 
{

    amrex::Print() << "\n##### Transport Solver #####\n\n";

    amrex::ParmParse pp_transport("transport");
    num_NS = 0;

    amrex::Vector<std::string> temp_vec;
    bool NS_pecified = pp_transport.queryarr("NS_names", temp_vec);

    int c=0;
    for (auto it: temp_vec)
    {
        if (std::find(vec_NS_names.begin(), vec_NS_names.end(), it) == vec_NS_names.end()) 
        {
           vec_NS_names.push_back(it);
           ++c;
        }
    }
    temp_vec.clear();

    amrex::Print() << "##### transport.NS_names: ";
    for (auto it: vec_NS_names) amrex::Print() << it << "  ";
    amrex::Print() << "\n";

}


void 
c_TransportSolver::InitData() 
{
    amrex::Print() << "\n##### TRANSPORT PROPERTIES #####\n\n";

    amrex::ParmParse pp_transport("transport");
    
    use_selfconsistent_potential = 0;
    use_negf = 0;
    pp_transport.query("use_selfconsistent_potential", use_selfconsistent_potential);
    pp_transport.query("use_negf", use_negf);
    amrex::Print() << "##### transport.use_selfconsistent_potential: " 
                   << use_selfconsistent_potential << "\n";
    amrex::Print() << "##### transport.use_negf: " 
                   << use_negf << "\n";

    std::string NS_gather_field_str = "phi";
    std::string NS_deposit_field_str = "charge_density";
    amrex::Real NS_initial_deposit_value = 0.;

    if(use_selfconsistent_potential) 
    {
        auto& rCode = c_Code::GetInstance();
        auto& rGprop = rCode.get_GeometryProperties();

        _geom = &rGprop.geom;
        _dm   = &rGprop.dm;
        _ba   = &rGprop.ba;

        pp_transport.get("NS_gather_field", NS_gather_field_str);
        amrex::Print() << "##### transport.NS_gather_field: " << NS_gather_field_str << "\n";
        if ( Evaluate_TypeOf_MacroStr(NS_gather_field_str) != 0 )
        {
            /*mf not defined in MacroscopicProperties*/
            amrex::Abort("NS_gather_field " + NS_gather_field_str + " not defined in Mprop.");
        } 

        pp_transport.get("NS_deposit_field", NS_deposit_field_str);
        amrex::Print() << "##### transport.NS_deposit_field: " << NS_deposit_field_str << "\n";
        if ( Evaluate_TypeOf_MacroStr(NS_deposit_field_str) != 0 )
        {
            /*mf not defined in MacroscopicProperties*/
            amrex::Abort("NS_deposit_field " + NS_deposit_field_str + " not defined in Mprop.");
        } 

        auto dxi = _geom->InvCellSizeArray();
        amrex::Real inv_vol = AMREX_D_TERM(dxi[0], *dxi[1], *dxi[2]);
        NS_initial_deposit_value = PhysConst::q_e*inv_vol;

        queryWithParser(pp_transport,"NS_initial_deposit_value", NS_initial_deposit_value);
        amrex::Print() << "##### transport.NS_initial_deposit_value: " << NS_initial_deposit_value << "\n";
    }

    std::string type;
    int c=0;
    for (auto name: vec_NS_names)
    {

        amrex::ParmParse pp_ns(name);
        pp_ns.get("type", type);
        amrex::Print() << "\n##### name: " << name << "\n";
        amrex::Print() << "##### type: " << type << "\n";

        switch (map_NSType_enum[type])
        {
            case s_NS::Type::CNT:
            {
                using T = c_CNT;
 
                vp_CNT.push_back(std::make_unique<c_Nanostructure<T>>(*_geom, *_dm, *_ba,
                                                                       name, 
                                                                       NS_gather_field_str, 
                                                                       NS_deposit_field_str, 
                                                                       NS_initial_deposit_value,
                                                                       use_selfconsistent_potential,
                                                                       use_negf
                                                                      ));
                break;
            }
            case s_NS::Type::Graphene:
            {
                using T = c_Graphene;
 
                vp_Graphene.push_back(std::make_unique<c_Nanostructure<T>>(*_geom, *_dm, *_ba,
                                                                            name, 
                                                                            NS_gather_field_str, 
                                                                            NS_deposit_field_str, 
                                                                            NS_initial_deposit_value,
                                                                            use_selfconsistent_potential,
                                                                            use_negf
                                                                           ));
                amrex::Abort("NS_type " + type + " is not yet defined.");
                break; 
            }
            case s_NS::Type::Silicon:
            {
                amrex::Abort("NS_type " + type + " is not yet defined.");
            }
            default:
            {
                amrex::Abort("NS_type " + type + " is not supported.");
            }
        }
        ++c;
    }

}

void 
c_TransportSolver::Solve() 
{

   for (int c=0; c < vp_CNT.size(); ++c)
   {
       if(use_selfconsistent_potential) 
       {
          vp_CNT[c]->GatherFromMesh();
          vp_CNT[c]->AverageFieldGatheredFromMesh();
          vp_CNT[c]->Write_AveragedGatherField();
       }

       //if(use_negf) 
       //{
       //    vp_CNT[c]->NEGF_Solve();
       //} 

   }

}