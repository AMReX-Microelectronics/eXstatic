#ifndef CNT_H_
#define CNT_H_

#include "Common_Properties.H"

#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>

#include <AMReX_Array.H>

#include <string>
#include <array>

class c_CNT
           :public c_Common_Properties
{

public:

    std::array<int,2> Get_Ring_Atom_IDs(int cell_id) 
    {
         std::array<int,2> id;
         id[0] = int((cell_id-1)/type_id[0]); //ring_id
         id[1] = (cell_id-1) % type_id[0]; //cell_id
         return id;  
    };

    const int get_1Dlayer_id(int global_par_id) 
    {
         return static_cast<int>(amrex::Math::floor((global_par_id-1)/type_id[0]));
    };

    const int get_atom_in_1Dlayer_id(int global_par_id) 
    {
         return (global_par_id-1) % type_id[0];
    };

    const int get_num_layers()
    {
         return rings_per_unitcell*num_unitcells;  
    };

    const int get_num_atoms_per_layer() 
    {
         return atoms_per_ring;
    };

    void Print_IDs(int cell_id) 
    {
         std::array<int,2> id;
         id = Get_Ring_Atom_IDs(cell_id);
         amrex::Print() << "ring id: " << id[0] << 
                           ", atom id: " << id[1] << "\n";

    };

    amrex::Array<int,2> type_id; //type_id[0]: m, type_id[1]: n for (m,n) tube
    amrex::Real acc;
    int rings_per_unitcell;
    int atoms_per_ring;
    amrex::Real gamma;

    virtual void ReadNanostructureProperties () override;

};

#endif
