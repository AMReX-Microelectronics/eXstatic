#include "BoundaryConditions.H"

#include "../../Utils/SelectWarpXUtils/WarpXUtil.H"
#include "../../Utils/SelectWarpXUtils/TextMsg.H"
#include "../../Utils/CodeUtils/CodeUtil.H"
#include "Code.H"
#include "GeometryProperties.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <ctype.h>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

using namespace amrex;

c_BoundaryConditions::c_BoundaryConditions ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_BoundaryConditions Constructor************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    DefineBoundaryTypeMap();
    some_robin_boundaries = false;
    ReadData();

#ifdef AMREX_USE_EB
        auto& rCode = c_Code::GetInstance();
        auto& rGprop = rCode.get_GeometryProperties();

        if(some_robin_boundaries == true and rGprop.embedded_boundary_flag==1)                      
        {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(some_robin_boundaries != false and  rGprop.embedded_boundary_flag !=0,
                                         "Use of robin boundaries with embedded boundaries is not supported/well tested at the moment. Consider changing those boundaries to Neumann.");
        }
#endif


#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_BoundaryConditions Constructor************************\n";
#endif
} 


c_BoundaryConditions::~c_BoundaryConditions ()
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t{************************c_BoundaryConditions Destructor()************************\n";
    amrex::Print() << "\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t}************************c_BoundaryConditions Destructor()************************\n";
#endif
}


void 
c_BoundaryConditions::ReadData()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t{************************c_BoundaryConditions::ReadData()************************\n";
    amrex::Print() << "\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    ReadBoundaryConditionsType();

    DefineMacroVariableVectorSizes();

    for (auto it: map_function_parser_name)
    {
        ReadBoundaryConditionsParser(it.first, it.second);
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t}************************c_BoundaryConditions::ReadData()************************\n";
#endif
}


void
c_BoundaryConditions::SortBoundaryTypeArrayString(const amrex::Vector<std::string>& bc_str, std::array< std::string, AMREX_SPACEDIM >& bcType, std::array< std::any, AMREX_SPACEDIM >& bcAny,  std::map<int,std::string>& map_bcAny)
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t\t{************************c_BoundaryConditions::SortBoundaryTypeArrayString(*)************************\n";
    amrex::Print() << "\t\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t\t\t\t";
#endif

    int c=0;
    for (auto str: bc_str)
    { 

        std::string first_three_letters = str.substr(0,3);
        bcType[c] = first_three_letters;
#ifdef PRINT_HIGH
        amrex::Print() << "\n" << prt << "string: " << str  << " first_three_letters: " << first_three_letters<< "\n";
#endif

        if(bcType[c] == "dir" or bcType[c] == "neu" or bcType[c] == "rob")
        { 
            std::string fourth_char = str.substr(3,1);
            std::string last_char = str.substr(str.length()-1);
#ifdef PRINT_HIGH
            amrex::Print() << prt << "fourth_char: " << fourth_char << "\n";
            amrex::Print() << prt << "last_char: " << last_char << "\n";
#endif

            if(fourth_char == "(" or last_char == ")" )
            {
                if(fourth_char == "(" and last_char == ")" )
                {
                }
                if(fourth_char == "(" and last_char != ")") //throw warning
                {
                    std::stringstream warnMsg;
                    warnMsg << "In the input file, specification of '" << str  << "' is missing a closed bracket\")\".\n"
                    << "Interpreting it as \n" << str.insert(str.length(),")") << "\n";
                    c_Code::GetInstance().RecordWarning("Boundary Conditions", warnMsg.str());

                }
                else if(fourth_char != "(" and last_char == ")") //throw warning
                {
                    std::stringstream warnMsg;
                    warnMsg << "In the input file, specification of '" << str  << "' is missing an open bracket\"(\".\n"
                    << "Interpreting it as \n" << str.insert(3,"(")  << "\n";
                    c_Code::GetInstance().RecordWarning("Boundary Conditions", warnMsg.str());

                }
                std::string bracketed_str = str.substr(4,str.length()-5);
               
#ifdef PRINT_HIGH
                amrex::Print() << prt << "bracketed_str: " << bracketed_str << "\n";
#endif
                std::string stripped_bracketed_str;
                if(bracketed_str.substr(0,1) == "-" or bracketed_str.substr(0,1) == "+") {
                   stripped_bracketed_str = bracketed_str.substr(1);
                }
                else {
                   stripped_bracketed_str = bracketed_str;
                }
                if(std::isdigit( *stripped_bracketed_str.c_str()) ) 
                {   
                    map_bcAny[c] = "inhomogeneous_constant"; 
                    bcAny[c] = std::stod(bracketed_str);
#ifdef PRINT_HIGH
                    amrex::Print() << prt << "inhomo constant: " << bracketed_str << "\n";
#endif
                }
                else 
                {
                    map_bcAny[c] = "inhomogeneous_function"; 
                    bcAny[c] = bracketed_str;
#ifdef PRINT_HIGH
                    amrex::Print() << prt << "inhomo function with parameter name: " << bracketed_str << "\n";
#endif
                }
            }
            else if(fourth_char != "(" and last_char != ")") 
            {
                map_bcAny[c] = "homogeneous"; 
                bcAny[c] = 0.0;
#ifdef PRINT_HIGH
                amrex::Print() << prt << "homo constant 0.0 " << "\n";
#endif
            }
        }
        else if(bcType[c] == "per") {
            map_bcAny[c] = "periodic"; 
        }
        else if(bcType[c] == "ref") {
            map_bcAny[c] = "reflect";
        }
        ++c;
    }

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t\t}************************c_BoundaryConditions::SortBoundaryTypeArrayString(*)************************\n";
#endif
}


void 
c_BoundaryConditions::ReadBoundaryConditionsType()
{ 

#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_BoundaryConditions::ReadBoundaryConditionsType()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t\t\t";
#endif

    amrex::Vector<amrex::Vector<std::string>> bc_str_2d(2);

    amrex::ParmParse pp_boundary("boundary");
    pp_boundary.queryarr("lo", bc_str_2d[0]);
    pp_boundary.queryarr("hi", bc_str_2d[1]);

    amrex::Print() << "\n##### BOUNDARY CONDITIONS #####\n\n";
    for (auto& i: bc_str_2d)
    {
        amrex::Print()  << "##### ";
        for (auto& j: i) 
        {
            amrex::Print()  << j << "  ";
        }
        amrex::Print() << "\n";
    }
        
    for (std::size_t i = 0; i < 2; ++i) 
    {
        SortBoundaryTypeArrayString(bc_str_2d[i], bcType_2d[i], bcAny_2d[i], map_bcAny_2d[i]);
    }

    bc_str_2d.clear();

    /*determine if there is a robin boundary*/
    for (std::size_t i = 0; i < 2; ++i) 
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j) 
        {
            if(bcType_2d[i][j] == "rob") {
                some_robin_boundaries = true;
            }  
        }
    }

    /* Make both boundaries periodic based on is_periodic */
    auto& rCode = c_Code::GetInstance();
    auto& rGprop = rCode.get_GeometryProperties();
    auto& is_periodic = rGprop.is_periodic;

    
    for (std::size_t idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
        if(is_periodic[idim] == 1 and  (map_bcAny_2d[0][idim] != "periodic" or  map_bcAny_2d[0][idim] != "periodic") ) 
        {
            std::stringstream warnMsg;
            warnMsg << "Note that domain.is_periodic is set to 1 (true) for direction "<< idim << " !\n"
                << "Therefore, the value set by boundary.lo/hi is ignored and both sides are assumed to be periodic. \n";
            c_Code::GetInstance().RecordWarning("Boundary Conditions", warnMsg.str());

            map_bcAny_2d[0][idim] = "periodic";    map_bcAny_2d[1][idim] = "periodic";
               bcType_2d[0][idim] = "per";            bcType_2d[1][idim] = "per";
                bcAny_2d[0][idim] = 0x0;               bcAny_2d[1][idim] = 0x0;

        }
    }

    /*Conversely, ensure a direction in is_periodic is set to be periodic if boundary.lo/hi are set to periodic.*/
    for (std::size_t idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
        bool is_periodic_flag = true;
        if ( ( bcType_2d[0][idim] == "per" and  bcType_2d[1][idim] != "per") or 
             ( bcType_2d[0][idim] != "per" and  bcType_2d[1][idim] == "per") or
             ( bcType_2d[0][idim] == "per" and  bcType_2d[1][idim] == "per" and is_periodic[idim] != 1) ) 
        {
             is_periodic_flag = false;
        }
         std::string warnMsg;
         warnMsg = "Both sides must be periodic, and 'domain.is_periodic', must be 1 for the periodic directions\
                    user intends to set. Current setup violates this for direction: " + std::to_string(idim);
         WARPX_ALWAYS_ASSERT_WITH_MESSAGE(is_periodic_flag == true, warnMsg);
    }


#ifdef PRINT_LOW
    for (std::size_t i = 0; i < 2; ++i) 
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j) 
        {
            
            amrex::Print()   << "\n" << prt << " boundary type, bcType_2d: " << bcType_2d[i][j] << "\n";
            amrex::Print()   << prt << " boundary map, map_bcAny_2d: " << map_bcAny_2d[i][j] << "\n";
            amrex::Print()   << prt << " (bcAny_2d) bracket ";

            process_std_any(bcAny_2d[i][j]);

            amrex::Print() << "\n";
        }
    }
    amrex::Print() << "\n";
#endif


    /* loop over map_bcAny_2d and fill in the map of function parser names and set number of function parser names. */
    int c=0;
    for (std::size_t i = 0; i < 2; ++i) 
    {
        for (std::size_t j = 0; j < AMREX_SPACEDIM; ++j) 
        {
            if(map_bcAny_2d[i][j] == "inhomogeneous_function") 
            {
                if( bcType_2d[i][j] == "rob") 
                {
                    std::string main_str = std::any_cast<std::string>(bcAny_2d[i][j]);

                    for(auto& str: robin_substrings) 
                    {
                        std::string str_with_subscript = main_str + str;
                        map_function_parser_name[ str_with_subscript ] = c; 
                        ++c;
                    }
                } 
                else 
                {
                    map_function_parser_name[ std::any_cast<std::string>(bcAny_2d[i][j]) ] = c; 
                    ++c;
                }
            }
        }
    }
    num_function_parsers = map_function_parser_name.size();

#ifdef PRINT_LOW
    amrex::Print() << "\n" << prt << "map_function_parser_name: \n";
    for (auto& i: map_function_parser_name)
    {
        amrex::Print() << prt << i.first << "  " << i.second << "\n";
    }
    amrex::Print() << "\n" << prt << "num_function_parsers: " << num_function_parsers << "\n";
    amrex::Print() << "\n\n";
#endif    

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_BoundaryConditions::ReadBoundaryConditionsType()************************\n\n";
#endif    
}


void 
c_BoundaryConditions::DefineMacroVariableVectorSizes()
{ 
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_BoundaryConditions::DefineMacroVariableVectorSizes()************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
#endif

    m_macro_str_function.resize(num_function_parsers);
    m_p_macro_parser.resize(num_function_parsers);

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_BoundaryConditions::DefineMacroVariableVectorSizes()************************\n";
#endif
}


void 
c_BoundaryConditions::ReadBoundaryConditionsParser(std::string macro_str, int macro_num)
{
#ifdef PRINT_NAME
    amrex::Print() << "\n\n\t\t\t\t\t{************************c_BoundaryConditions::ReadBoundaryConditionsParser(*)************************\n";
    amrex::Print() << "\t\t\t\t\tin file: " << __FILE__ << " at line: " << __LINE__ << "\n";
    std::string prt = "\t\t\t\t\t";
#endif


    ParmParse pp_boundary("boundary");

    //std::string function_str;
    //#ifdef TIME_DEPENDENT
    //    function_str = "_function(x,y,z,t)";
    //#else
    //    function_str = "_function(x,y,z)";
    //#endif

    std::string macro_functionXYZ = macro_str + "_function";
    bool specified = false;

    if (pp_boundary.query( macro_functionXYZ.c_str(), m_macro_str_function[macro_num]) ) {
        specified = true;
    }
    if(!specified) {
        std::string warnMsg = "Boundary Conditions: function parser '" + macro_functionXYZ + "' is not specified in the input file.\n";
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(specified==true, warnMsg);
    } 
    else 
    {
        Store_parserString(pp_boundary, macro_functionXYZ.c_str(),  m_macro_str_function[macro_num]);

        m_p_macro_parser[macro_num] = std::make_unique<amrex::Parser>(
                                           makeParser( m_macro_str_function[macro_num], m_parser_varname_vector));
    }

#ifdef PRINT_LOW
    if(specified) {
       amrex::Print() << "\n" << prt << "Reading macro: " << macro_str << " with number: " << macro_num << "\n";
       amrex::Print() << prt << "function_parser_name: " << macro_functionXYZ << "\n";
       amrex::Print() << prt << "function_parser: " <<  m_macro_str_function[macro_num] << "\n\n";
    }
#endif

#ifdef PRINT_NAME
    amrex::Print() << "\t\t\t\t\t}************************c_BoundaryConditions::ReadBoundaryConditionsParser(*)************************\n\n";
#endif    
}
