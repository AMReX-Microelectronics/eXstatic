#ifndef CNT_H_
#define CNT_H_

#include "NEGF_Common.H"
#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>

#include <string>
#include <array>


class c_CNT
           :public c_Common_Properties<ComplexType[NUM_MODES]>
{

public:
    using BlkType = ComplexType[NUM_MODES];

    std::array<int,2> Get_Ring_Atom_IDs(int cell_id) 
    {
         std::array<int,2> id;
         id[0] = int((cell_id-1)/type_id[0]); //ring_id
         id[1] = (cell_id-1) % type_id[0]; //cell_id
         return id;  
    };

    //AMREX_GPU_HOST_DEVICE static const int get_1Dlayer_id(int global_par_id) 
    //{
    //     return static_cast<int>(amrex::Math::floor((global_par_id-1)/type_id[0]));
    //};

    struct Get1DLayerID {
       AMREX_GPU_HOST_DEVICE int operator() (int global_par_id) const 
       { 
	   return static_cast<int>(amrex::Math::floor((global_par_id-1)/type_id[0])); 
       }
       amrex::Array<int,2> type_id;
    };

    static Get1DLayerID get_1Dlayer_id() {
         return Get1DLayerID{type_id}; 
    }

    struct GetAtomIn1DLayerID {
       AMREX_GPU_HOST_DEVICE int operator() (int global_par_id) const 
       { 
	   return static_cast<int>((global_par_id-1) % type_id[0]); 
       }
       amrex::Array<int,2> type_id;
    };

    static GetAtomIn1DLayerID get_atom_in_1Dlayer_id() {
         return GetAtomIn1DLayerID{type_id}; 
    }

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


    static amrex::Array<int,2> type_id; //type_id[0]: m, type_id[1]: n for (m,n) tube
    amrex::Real acc;
    int rings_per_unitcell;
    int atoms_per_ring;
    amrex::Real gamma;
    amrex::Array<const int, NUM_MODES> mode_arr = {6,11};
    amrex::Array<const int, NUM_MODES> mode_degen_arr = {2,2};
    MatrixBlock<BlkType> beta;

    virtual int get_Hsize() override { return rings_per_unitcell*num_unitcells; }
    virtual int get_offDiag_repeatBlkSize() override { return 2; }

    virtual void ReadNanostructureProperties () override;
    virtual void Define_MPI_BlkType () override;
    virtual void AllocateArrays() override;
    virtual void ConstructHamiltonian () override;
    virtual void Define_EnergyLimits () override;
    virtual void Define_IntegrationPaths () override;

    virtual void Define_ContactInfo() override;
    virtual void Compute_SurfaceGreensFunction (MatrixBlock<BlkType>& gr, const ComplexType EmU) override;

    virtual void Define_SelfEnergy() override;
    //virtual void DeallocateArrays() override;
    //virtual void ComputeChargeDensity() override;


    AMREX_GPU_HOST_DEVICE
    ComplexType get_beta(int J);

    AMREX_GPU_HOST_DEVICE
    void get_Sigma(MatrixBlock<BlkType>& Sigma, const ComplexType EmU);
    //void get_Sigma(MatrixBlockCNT& Sigma, const ComplexType EmU);

};

#endif