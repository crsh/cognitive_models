gfortran Number_generators.f90 MinervaAL_tools.f90 Acquisition_1.f90 -o acquisition_1
gfortran Number_generators.f90 MinervaAL_tools.f90 Acquisition_067.f90 -o acquisition_067
gfortran Number_generators.f90 MinervaAL_tools.f90 Acquisition_033.f90 -o acquisition_033

gfortran Number_generators.f90 MinervaAL_tools.f90 Reacquisition_1.f90 -o reacquisition_1
gfortran Number_generators.f90 MinervaAL_tools.f90 Reacquisition_067.f90 -o reacquisition_067
gfortran Number_generators.f90 MinervaAL_tools.f90 Reacquisition_033.f90 -o reacquisition_033

gfortran Number_generators.f90 MinervaAL_tools.f90 Reacquisition_control_1.f90 -o reacquisition_control_1
gfortran Number_generators.f90 MinervaAL_tools.f90 Reacquisition_control_067.f90 -o reacquisition_control_067
gfortran Number_generators.f90 MinervaAL_tools.f90 Reacquisition_control_033.f90 -o reacquisition_control_033

./acquisition_1
./acquisition_067
./acquisition_033

./reacquisition_1
./reacquisition_067
./reacquisition_033

./reacquisition_control_1
./reacquisition_control_067
./reacquisition_control_033