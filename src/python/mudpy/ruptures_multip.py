"""
The function generate_ruptures can be used in place of
fakequakes.generate_ruptures to generate a list of ruptures in parallel
using the multiprocessing module rather than mpi4py.

This allows setting the random seed for each rupture separately in
order to make the code reproducible.

Requires the class fakequakes_utils.Realization to specify input
data for each desired realization.
"""

def generate_ruptures(home,project_name,run_name,fault_name,slab_name,mesh_name,
		load_distances,distances_name,UTM_zone,realizations,model_name,hurst,Ldip,
		Lstrike,num_modes,rake,rise_time_depths,time_epi,
		max_slip,source_time_function,lognormal,slip_standard_deviation,scaling_law,ncpus,
		force_magnitude=False,force_area=False,mean_slip_name=None,hypocenter=None,
		slip_tol=1e-2,force_hypocenter=False,no_random=False,shypo=None,use_hypo_fraction=True,
		shear_wave_fraction_shallow=0.49,shear_wave_fraction_deep=0.8,max_slip_rule=True):


    from numpy import load,save,genfromtxt,log10,cos,sin,deg2rad,savetxt,zeros,where
    from time import gmtime, strftime
    from numpy.random import shuffle
    from numpy import random
    from mudpy import fakequakes
    from obspy import UTCDateTime
    from obspy.taup import TauPyModel
    import geopy.distance
    import warnings
    from numpy import ceil
    from os import environ
    import numpy as np
    import multiprocessing
    from multiprocessing import Process, current_process
    import time


    # Need to make tauPy file  (Could move this to initialization?)
    print('Creating tauPy file...')
    vel_mod_file=home+project_name+'/structure/'+model_name
    fakequakes.build_TauPyModel(home,project_name,vel_mod_file,
                                background_model='PREM')


    
    print("multiprocessing: Using %s CPUs" % ncpus)

    #I don't condone it but this cleans up the warnings
    warnings.filterwarnings("ignore")

            
    if time_epi=='None':
        time_epi=None
    else:
        time_epi=UTCDateTime(time_epi)
    #rise_time_depths=[rise_time_depths0,rise_time_depths1]

    

    #Should I calculate or load the distances?
    if load_distances==1:  
        Dstrike=load(home+project_name+'/data/distances/'+distances_name+'.strike.npy')
        Ddip=load(home+project_name+'/data/distances/'+distances_name+'.dip.npy')
    else:
        Dstrike,Ddip=fakequakes.subfault_distances_3D(home,project_name,fault_name,slab_name,UTM_zone)
        save(home+project_name+'/data/distances/'+distances_name+'.strike.npy',Dstrike)
        save(home+project_name+'/data/distances/'+distances_name+'.dip.npy',Ddip)
    

    #Read fault and prepare output variable
    whole_fault=genfromtxt(home+project_name+'/data/model_info/'+fault_name)
    
    #Get structure model
    vel_mod_file=home+project_name+'/structure/'+model_name
    
    #Get TauPyModel
    velmod = TauPyModel(model=home+project_name+'/structure/'+model_name.split('.')[0])

    #Now loop over the number of realizations
    #realization=0  # not used, nor is Nrealizations
    
    
    def generate_one_rupture(realization):
        current_target_Mw = realization.target_Mw
        rname = realization.rname
        rseed = realization.rseed
        if rseed is None:
            # generate random seed if missing:
            rseed = random.randint(100000,10000000,1)
            realization.rseed = rseed
            print('*** Generated random seed, but may not be distinct')
        random.seed(rseed)  # initialize random number generator
        
        #Prepare output
        fault_out=zeros((len(whole_fault),14))
        fault_out[:,0:8]=whole_fault[:,0:8]
        fault_out[:,10:12]=whole_fault[:,8:]   
        
        #Sucess criterion
        success=False
        while success==False:
            #Select only a subset of the faults based on magnitude scaling
            
            ifaults,hypo_fault,Lmax,Wmax,Leff,Weff=fakequakes.select_faults(whole_fault,Dstrike,Ddip,current_target_Mw,num_modes,scaling_law,
                                force_area,no_shallow_epi=False,no_random=no_random,subfault_hypocenter=shypo,use_hypo_fraction=use_hypo_fraction)
            fault_array=whole_fault[ifaults,:]
            Dstrike_selected=Dstrike[ifaults,:][:,ifaults]
            Ddip_selected=Ddip[ifaults,:][:,ifaults]
            
            #Determine correlation lengths from effective length.width Leff and Weff
            if Lstrike=='MB2002': #Use scaling
                #Ls=10**(-2.43+0.49*target_Mw)
                Ls=2.0+(1./3)*Leff
            elif Lstrike=='auto':
                Ls=17.7+0.34*Leff
            else:
                Ls=Lstrike
            if Ldip=='MB2002': #Use scaling
                #Ld=10**(-1.79+0.38*target_Mw)
                Ld=1.0+(1./3)*Weff
            elif Ldip=='auto':
                Ld=6.8+0.4*Weff
            else:
                Ld=Ldip
            
            #Get the mean uniform slip for the target magnitude
            if mean_slip_name==None:
                mean_slip,mu=fakequakes.get_mean_slip(current_target_Mw,fault_array,vel_mod_file)
            else:
                foo,mu=fakequakes.get_mean_slip(current_target_Mw,fault_array,vel_mod_file)
                mean_fault=genfromtxt(mean_slip_name)
                mean_slip=(mean_fault[:,8]**2+mean_fault[:,9]**2)**0.5
                
                #keep onlt faults that have man slip inside the fault_array seelcted faults
                mean_slip=mean_slip[ifaults]
                
                #get the area in those selected faults
                area=fault_array[:,-2]*fault_array[:,-1]
                
                #get the moment in those selected faults
                moment_on_selected=(area*mu*mean_slip).sum()
                
                #target moment
                target_moment=10**(1.5*current_target_Mw+9.1)
                
                #How much do I need to upscale?
                scale_factor=target_moment/moment_on_selected
                
                #rescale the slip
                mean_slip = mean_slip*scale_factor
                
                
                #Make sure mean_slip has no zero slip faults
                izero=where(mean_slip==0)[0]
                mean_slip[izero]=slip_tol
            
            #Get correlation matrix
            C=fakequakes.vonKarman_correlation(Dstrike_selected,Ddip_selected,Ls,Ld,hurst)
            
            # Lognormal or not?
            if lognormal==False:
                #Get covariance matrix
                C_nonlog=fakequakes.get_covariance(mean_slip,C,current_target_Mw,fault_array,vel_mod_file,slip_standard_deviation) 
                #Get eigen values and eigenvectors
                eigenvals,V=fakequakes.get_eigen(C_nonlog)
                #Generate fake slip pattern
                rejected=True
                while rejected==True:
#                        slip_unrectified,success=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip,max_slip,lognormal=False,seed=kfault)
                    slip_unrectified,success=fakequakes.make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip,max_slip,lognormal=False,seed=12345)
                    slip,rejected,percent_negative=fakequakes.rectify_slip(slip_unrectified,percent_reject=13)
                    if rejected==True:
                        print('... ... ... negative slip threshold exceeeded with %d%% negative slip. Recomputing...' % (percent_negative))
            else:
                #Get lognormal values
                C_log,mean_slip_log=fakequakes.get_lognormal(mean_slip,C,current_target_Mw,fault_array,vel_mod_file,slip_standard_deviation)               
                #Get eigen values and eigenvectors
                eigenvals,V=fakequakes.get_eigen(C_log)
                #Generate fake slip pattern
#                    slip,success=make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip_log,max_slip,lognormal=True,seed=kfault)
                slip,success=fakequakes.make_KL_slip(fault_array,num_modes,eigenvals,V,mean_slip_log,max_slip,lognormal=True,seed=12345)
        
            #Slip pattern sucessfully made, moving on.
            #Rigidities
            foo,mu=fakequakes.get_mean_slip(current_target_Mw,whole_fault,vel_mod_file)
            fault_out[:,13]=mu
            
            #Calculate moment and magnitude of fake slip pattern
            M0=sum(slip*fault_out[ifaults,10]*fault_out[ifaults,11]*mu[ifaults])
            Mw=(2./3)*(log10(M0)-9.1)
            
            #Check max_slip_rule
            if max_slip_rule==True:
                
                max_slip_from_rule=10**(-4.94+0.71*Mw) #From Allen & Hayes, 2017
                max_slip_tolerance = 3
                
                if slip.max() > max_slip_tolerance*max_slip_from_rule:
                    success = False
                    print('... ... ... max slip condition violated max_slip_rule, recalculating...')
            
            #Force to target magnitude
            if force_magnitude==True:
                M0_target=10**(1.5*current_target_Mw+9.1)
                M0_ratio=M0_target/M0
                #Multiply slip by ratio
                slip=slip*M0_ratio
                #Recalculate
                M0=sum(slip*fault_out[ifaults,10]*fault_out[ifaults,11]*mu[ifaults])
                Mw=(2./3)*(log10(M0)-9.1)
                
            #check max_slip again
            if slip.max() > max_slip:
                success=False
                print('... ... ... max slip condition violated due to force_magnitude=True, recalculating...')
        
        
        #Get stochastic rake vector
        stoc_rake=fakequakes.get_stochastic_rake(rake,len(slip))
        
        #Place slip values in output variable
        fault_out[ifaults,8]=slip*cos(deg2rad(stoc_rake))
        fault_out[ifaults,9]=slip*sin(deg2rad(stoc_rake))
        
        #Move hypocenter to somewhere with a susbtantial fraction of peak slip
#            slip_fraction=0.25
#            islip=where(slip>slip.max()*slip_fraction)[0]
#            shuffle(islip) #randomize
#            hypo_fault=ifaults[islip[0]] #select first from randomized vector
        
        #Calculate and scale rise times
        rise_times=fakequakes.get_rise_times(M0,slip,fault_array,rise_time_depths,stoc_rake)
        
        #Place rise_times in output variable
        fault_out[:,7]=0
        fault_out[ifaults,7]=rise_times
        
        #Calculate rupture onset times
        if force_hypocenter==False: #Use random hypo, otehrwise force hypo to user specified
            hypocenter=whole_fault[hypo_fault,1:4]
        else:
            hypocenter=whole_fault[shypo,1:4]
            
        t_onset=fakequakes.get_rupture_onset(home,project_name,slip,fault_array,model_name,hypocenter,rise_time_depths,
                                             M0,velmod,shear_wave_fraction_shallow=shear_wave_fraction_shallow,
                                             shear_wave_fraction_deep=shear_wave_fraction_deep)
        fault_out[:,12]=0
        fault_out[ifaults,12]=t_onset
        
        #Calculate location of moment centroid
        centroid_lon,centroid_lat,centroid_z=fakequakes.get_centroid(fault_out)
        
        #Calculate average risetime
        rise = fault_out[:,7]
        avg_rise = np.mean(rise)
        
        # Calculate average rupture velocity
        lon_array = fault_out[:,1]
        lat_array = fault_out[:,2]
        vrupt = []
        
        for i in range(len(fault_array)):
            if t_onset[i] > 0:
                r = geopy.distance.geodesic((hypocenter[1], hypocenter[0]), (lat_array[i], lon_array[i])).km
                vrupt.append(r/t_onset[i])
        
        avg_vrupt = np.mean(vrupt)
        
        #Write to file
    
        outfile=home+project_name+'/output/ruptures/'+rname+'.rupt'
        savetxt(outfile,fault_out,fmt='%d\t%10.6f\t%10.6f\t%8.4f\t%7.2f\t%7.2f\t%4.1f\t%5.2f\t%5.2f\t%5.2f\t%10.2f\t%10.2f\t%5.2f\t%.6e',header='No,lon,lat,z(km),strike,dip,rise,dura,ss-slip(m),ds-slip(m),ss_len(m),ds_len(m),rupt_time(s),rigidity(Pa)')
        
        #Write log file
        logfile=home+project_name+'/output/ruptures/'+rname+'.log'
        f=open(logfile,'w')
        f.write('Scenario calculated at '+strftime("%Y-%m-%d %H:%M:%S", gmtime())+' GMT\n')
        f.write('Project name: '+project_name+'\n')
        f.write('Run name: '+run_name+'\n')
        f.write('Realization name: '+rname+'\n')
        f.write('Random seed: '+str(rseed)+'\n')
        f.write('Velocity model: '+model_name+'\n')
        f.write('No. of KL modes: '+str(num_modes)+'\n')
        f.write('Hurst exponent: '+str(hurst)+'\n')
        f.write('Corr. length used Lstrike: %.2f km\n' % Ls)
        f.write('Corr. length used Ldip: %.2f km\n' % Ld)
        f.write('Slip std. dev.: %.3f km\n' % slip_standard_deviation)
        f.write('Maximum length Lmax: %.2f km\n' % Lmax)
        f.write('Maximum width Wmax: %.2f km\n' % Wmax)
        f.write('Effective length Leff: %.2f km\n' % Leff)
        f.write('Effective width Weff: %.2f km\n' % Weff)
        f.write('Target magnitude: Mw %.4f\n' % current_target_Mw)
        f.write('Actual magnitude: Mw %.4f\n' % Mw)
        f.write('Hypocenter (lon,lat,z[km]): (%.6f,%.6f,%.2f)\n' %(hypocenter[0],hypocenter[1],hypocenter[2]))
        f.write('Hypocenter time: %s\n' % time_epi)
        f.write('Centroid (lon,lat,z[km]): (%.6f,%.6f,%.2f)\n' %(centroid_lon,centroid_lat,centroid_z))
        f.write('Source time function type: %s\n' % source_time_function)
        f.write('Average Risetime (s): %.2f\n' % avg_rise)
        f.write('Average Rupture Velocity (km/s): %.2f' % avg_vrupt)
        f.close()
        
        #realization+=1
        
    def generate_all_ruptures(realizations, nprocs):
        
        print("\n%s realizations will be generated on %s processors" \
                % (len(realizations),nprocs))
                
        if nprocs <= 1:
            # serial mode, no need to use multiprocess:
            
            for realization in realizations:
                print("Now generating this rupture: %s"  % realization.rname)
                generate_one_rupture(realization)
        else:
            
            # Split up cases between nprocs processors:
            # caseproc[i] will be a list of all cases to be handled by processor i
            caseproc = [[] for p in range(nprocs)]
            for i in range(len(realizations)):
                caseproc[i%nprocs].append(realizations[i])

            
            abort_time = 5
            print("You have %s seconds to abort..." % abort_time)

            time.sleep(abort_time) # give time to abort
            

            def run_cases(procnum):
                # loop over all cases assigned to one processor:
                p = current_process()
                time.sleep(2)
                message =  "\nProcess %s will run these cases:\n" % p.pid
                for realization in caseproc[procnum]:
                    message = message + "  %s \n" % realization.rname
                print(message) # constructed first to avoid interleaving prints
                time.sleep(2)

                for realization in caseproc[procnum]:
                    print("\nProcess %s now generating this rupture: %s, seed=%s" \
                            % (p.pid, realization.rname, realization.rseed))
                    generate_one_rupture(realization)
                    
            plist = [Process(target=run_cases, args=(p,)) for p in range(nprocs)]

            # Handle termination of subprocesses in case user terminates the main process
            def terminate_procs():
                for P in plist:
                    if P.is_alive():
                        P.terminate()
            # If this is registered, all processes will die if any are killed, e.g.
            # if ctrl-C received, but also if python exception encountered in
            # postprocessing of any run...
            #atexit.register(terminate_procs)

            print("\n=========================================\n")
            for P in plist:
                P.start()
                print("Starting process: ",P.pid)
            for P in plist:
                P.join()


    generate_all_ruptures(realizations, ncpus)

    # write names to ruptures.list (append to old list)
    fname = home+project_name+'/data/ruptures.list'
    with open(fname, 'a') as f:
        for realization in realizations:
            rname = realization.rname
            f.write('%s\n' % rname)
    print('Appended %i realizations to %s' % (len(realizations), fname))
    print('Ruptures are in %s' % home+project_name+'/ruptures')
        
