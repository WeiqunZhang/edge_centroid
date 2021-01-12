#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB_utils.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        int n_cell = 128;
        int max_grid_size = 64;

        // read parameters
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
        }

        Geometry geom;
        {
            RealBox rb({-1.0,-1.0,-1.0}, {1.0,1.0,1.0}); // physical domain
            Array<int,AMREX_SPACEDIM> is_periodic{false, false, false};
            Box domain(IntVect(0), IntVect(n_cell-1));
            geom.define(domain, rb, CoordSys::cartesian, is_periodic);
        }

        EB2::SphereIF sphere(0.5, {0.0,0.0,0.0}, false);
        auto gshop = EB2::makeShop(sphere);
        EB2::Build(gshop, geom, 0, 0);

        BoxArray ba(geom.Domain());
        ba.maxSize(max_grid_size);

        DistributionMapping dm{ba};

        std::unique_ptr<EBFArrayBoxFactory> eb_fact = makeEBFabFactory(geom,ba,dm,
                                                                       {2,2,2},
                                                                       EBSupport::full);
        Array<MultiFab,AMREX_SPACEDIM> edge_centroid;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            IntVect index_type{1};
            index_type[idim] = 0;  // All nodal except idim-direction
            edge_centroid[idim].define(amrex::convert(ba,index_type), dm, 1, 1, MFInfo{}, *eb_fact);
        }

        amrex::FillEdgeCentroid(GetArrOfPtrs(edge_centroid));

        for (MFIter mfi(edge_centroid[0]); mfi.isValid(); ++mfi) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const& edge = edge_centroid[idim].const_array(mfi);
                amrex::LoopOnCpu(Box(edge), [=] (int i, int j, int k)
                {
                    if (edge(i,j,k) == Real(-1.0)) {
                        // This edge is all covered
                    } else if (edge(i,j,k) == Real(1.0)) {
                        // This edge is all open
                    } else {
                        // This edge is cut.
                        Real edge_cent = edge(i,j,k); // edge centroid: (-0.5,0.5)
                        if (edge_cent < Real(0.0)) {
                            // The right side is covered.
                            Real cut_position = Real(2.0)*edge_cent + Real(0.5); // (-0.5,0.5)
                        } else {
                            // The left side is covered
                            Real cut_position = Real(2.0)*edge_cent - Real(0.5); // (-0.5,0.5)
                        }
                        Real cut_poisition = Real(2.0)*edge_cent
                            - std::copysign(Real(0.5), edge_cent); // (-0.5, 0.5)
                    }
                });
            }
        }
    }

    amrex::Finalize();
}
