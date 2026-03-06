cmake_minimum_required(VERSION 3.16)
project(MC_simple_stage1 C)

set(DSFMT_MEXP 19937)

add_executable(MC_simple
    main.c
    memory.c
    input_parser.c
    lattice.c
    mc_update.c
    exmc.c
    dSFMT.c
)
target_compile_definitions(MC_simple PRIVATE
    NO_MPI
    DSFMT_MEXP=${DSFMT_MEXP}
)
target_compile_options(MC_simple PRIVATE -O2)
target_link_libraries(MC_simple PRIVATE m)
