!------------------------------------
! Acquisition and Extinction 
!------------------------------------
program Acquisition_Extinction
use Number_generators
use MinervaAL_tools
implicit none

	integer, parameter		::	N_cues = 10,							&
								N_field = 20,							&
								N_features = N_cues*N_field,			&
								N_trials = 100,							&
								N_subjects = 100,						&								
								N_phases = 2
 
	real, parameter			::	L = 1

	real					::	Echo(N_features),						&
								Probe(N_features),						&	
								Cue_matrix(N_features, N_cues),			&
								Memory(N_features, N_trials*N_phases),	&
								Summary(N_trials*N_phases, N_subjects)
														
	integer					::	i, 										&
								j, 										& 
								k, 										& 
								v, 										& 
								N_traces


 	call RandSeed
	call Get_stimuli(Cue_matrix, N_cues)

	Summary = 0.0

	do i = 1, N_subjects
	N_traces = 0
	Memory = 0.0	  
		do j = 1, N_phases
		do k = 1, N_trials
					
			!---------------------------------------------------	
			! Construct the probe to be relevant for the 
			! current learning phase: 1 = A, 9 = Context, 
			! and 10 = Outcome (X). For example, j == 1 is an
			! A+ trial whereas j == 2 is an A- trial
			!---------------------------------------------------
			Probe = 0.0	
			if (j == 1) Probe(:) = Cue_matrix(:, 1) + Cue_matrix(:, 9) + Cue_matrix(:, 10) 	! Acquisition
			if (j == 2) Probe(:) = Cue_matrix(:, 1) + Cue_matrix(:, 9)						! Extinction
				
			!---------------------------------------------------	
			! Get the echo for the probe with 0.001 noise added
			!---------------------------------------------------
			call Get_echo(Echo, Probe, Memory, N_traces, 0.001)

			!---------------------------------------------------	
			! Increment N_traces, store the response strength,
			! and encode memory for the trial
			!---------------------------------------------------
			N_traces = N_traces + 1
			Summary(N_traces, i) = Sum(Echo(181:200))/N_field	! A shorthand of formula 5 in Jamieson et al. (2012)
			do v = 1, N_features
				if (flat(1.0) < L) Memory(v, N_traces) = Probe(v) - Echo(v)
			enddo

		enddo 
		enddo 
	enddo

	!---------------------------------------------------	
	! Write the results of the simulation to a file
	! as a matrix of N_subjects rows by N_trials columns
	!---------------------------------------------------
	Open(1, file='results/Acquisition_extinction_1.txt')
	write(1,*)
	write(1,'(A, F5.2)') 'Learning rate = ', L
	write(1,*)
	do i = 1, N_subjects
		write(1, '(A2, I4, 1000F8.4)') 'Ss ', i, Summary(:, i)
	enddo
	write(1,*)
	write(1,'(A6)',advance='no') 'M' 
	do i = 1, N_traces
		write(1, '(1000F8.4)',advance='no') Mean(Summary(i,:), N_subjects)
	enddo
	write(1,*)
	write(1,'(A6)',advance='no') 'SEM' 
	do i = 1, N_traces
		write(1, '(1000F8.4)',advance='no') SEM(Summary(i,:), N_subjects)
	enddo
	write(1,*)	
	Close(1)

END PROGRAM Acquisition_extinction
