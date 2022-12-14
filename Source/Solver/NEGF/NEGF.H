#ifndef NEGF_H_
#define NEGF_H_

#include "NEGF_fwd.H"
#include "Nanostructure.H"

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

#include <string>

struct s_NS {
    enum Type { CNT, Graphene, Silicon };
};

class 
c_NEGFSolver
{

public:
    amrex::Vector<std::unique_ptr<c_Nanostructure<c_CNT>>> vp_CNT;
    amrex::Vector<std::unique_ptr<c_Nanostructure<c_Graphene>>> vp_Graphene;
    amrex::Vector<std::unique_ptr<c_Nanostructure<c_Silicon>>> vp_Silicon;

    c_NEGFSolver();
    ~c_NEGFSolver();

    void InitData();
    void Solve();

private:
    int num_NS=0;
    const amrex::Geometry * _geom;
    const amrex::BoxArray * _ba;
    const amrex::DistributionMapping * _dm;
    //const amrex::GpuArray<int, AMREX_SPACEDIM> * _n_cell;
    amrex::Vector< std::string > vec_NS_names;

    void ReadData();

    std::map<std::string, s_NS::Type>
    map_NSType_enum = 
        {
           {"carbon nanotube" , s_NS::Type::CNT},
           {"Carbon Nanotube" , s_NS::Type::CNT},
           {"CNT"             , s_NS::Type::CNT},
           {"graphene"        , s_NS::Type::Graphene},
           {"Graphene"        , s_NS::Type::Graphene},
           {"silicon"         , s_NS::Type::Silicon},
           {"Silicon"         , s_NS::Type::Silicon}
        };

    std::string NS_gather_field_str;
    std::string NS_deposit_field_str;
    amrex::Real NS_initial_deposit_value;

};

#endif
