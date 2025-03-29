# axisymmetric_membranes_contact_force

g++ -I /usr/local/include/eigen3/ -I /opt/intel/oneapi/mkl/2025.0/include/ main.cpp world.cpp setInput.cpp timeStepper.cpp inertialForce.cpp externalGravityForce.cpp dampingForce.cpp elasticPlate.cpp externalPressureForce.cpp hyperElasticMM.cpp -lGL -lglut -lGLU -L /opt/intel/oneapi/mkl/2025.0/lib -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L /opt/intel/oneapi/compiler/2025.0/lib -liomp5 -llapack -lgfortran -fopenmp -lpthread -lm -Ofast -o simDER
