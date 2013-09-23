#============================================================================================
global CURRENT_VRSION
CURRENT_VRSION = "v.2.03"
import matplotlib
#matplotlib.use('WXAgg')

#from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas \

import sys,pylab,scipy,os
##try:
##    import pmag
##except:
##    pass
try:
    import thellier_gui_preferences
except:
    pass
import copy
import stat
import subprocess
import time
#import wx
#import wx.grid
import random
import copy
from pylab import *
from scipy.optimize import curve_fit
import wx.lib.agw.floatspin as FS
try:
    from mpl_toolkits.basemap import Basemap, shiftgrid
except:
    pass

#from matplotlib.backends.backend_wx import NavigationToolbar2Wx

import thellier_consistency_test

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
matplotlib.rc('axes', labelsize=8) 
matplotlib.rcParams['savefig.dpi'] = 300.

rcParams.update({"svg.embed_char_paths":False})
rcParams.update({"svg.fonttype":'none'})



#============================================================================================



    
class Arai_GUI():
    """ The main frame of the application
    """
    title = "PmagPy Thellier GUI %s"%CURRENT_VRSION
    
    def __init__(self, magic_file = "magic_measurements.txt"):
        print " calling __init__ Arai_gui instance"
        self.redo_specimens={}
        self.currentDirectory = "/Users/nebula/Python/SPD_project"
        self.WD = "/Users/nebula/Python/SPD_project"
        accept_new_parameters_default,accept_new_parameters_null=self.get_default_criteria()    # inialize Null selecting criteria
        self.accept_new_parameters_null=accept_new_parameters_null
        self.accept_new_parameters_default=accept_new_parameters_default
        #self.accept_new_parameters=copy.deepcopy(accept_new_parameters_default)
        preferences=[]
        self.dpi = 100
        self.magic_file="/Users/nebula/Python/SPD_project/" + magic_file
        self.preferences=preferences
        # inialize selecting criteria
        accept_new_parameters=self.read_criteria_from_file(self.WD+"/pmag_criteria.txt")          
        self.accept_new_parameters=accept_new_parameters
        #self.accept_new_parameters=accept_new_parameters
        self.Data,self.Data_hierarchy,self.Data_info={},{},{}
        self.MagIC_directories_list=[]

        self.Data_info=self.get_data_info() # get all ages, locations etc. (from er_ages, er_sites, er_locations)
        print "arai_GUI initialization calling self.get_data()"
        self.Data,self.Data_hierarchy=self.get_data() # Get data from magic_measurements and rmag_anistropy if exist.

        self.Data_samples={}
        self.last_saved_pars={}
        self.specimens=self.Data.keys()         # get list of specimens
        self.specimens.sort()                   # get list of specimens

        self.get_previous_interpretation() # get interpretations from pmag_specimens.txt
        print "data info: ", self.Data_info
    

    #----------------------------------------------------------------------
    #----------------------------------------------------------------------        



    #----------------------------------------------------------------------        


    def read_criteria_from_file(self,criteria_file):

        """
        Try to read pmag_criteria.txt from working directory
        return a full list of acceptance criteria
        If cant read file rerurn the default values 
        """
        print "calling read_criteria_from_file"
        # initialize Null and Default 
        default_acceptance_criteria,null_acceptance_criteria=self.get_default_criteria()        # Replace with new parametrs
        replace_acceptance_criteria={}
        for key in null_acceptance_criteria:
            replace_acceptance_criteria[key]=null_acceptance_criteria[key]
        try:
            fin=open(criteria_file,'rU')
            line=fin.readline()
            line=fin.readline()
            header=line.strip('\n').split('\t')
            for L in fin.readlines():
                line=L.strip('\n').split('\t')
                for i in range(len(header)):

                        if header[i] in self.high_threshold_velue_list + ['anisotropy_alt']:
                            try:
                                if float(line[i])<100:
                                    replace_acceptance_criteria[header[i]]=float(line[i])
                            except:
                                pass
                        if header[i] in self.low_threshold_velue_list:
                            try:
                                if float(line[i])>0.01:
                                    replace_acceptance_criteria[header[i]]=float(line[i])
                            except:
                                pass

                        # scat parametr (true/false)
                        if header[i] == 'specimen_scat' and ( line[i]=='True' or line[i]=='TRUE' or line[i]==True):
                                replace_acceptance_criteria['specimen_scat']=True
                        if header[i] == 'specimen_scat' and ( line[i]=='False' or line[i]=='FALSE' or line[i]==False):
                                replace_acceptance_criteria['specimen_scat']=False

                        # aniso parametr (true/false)
                        if header[i] == 'check_aniso_ftest' and ( line[i]=='True' or line[i]=='TRUE' or line[i]==True):
                                replace_acceptance_criteria['check_aniso_ftest']=True
                        if header[i] == 'check_aniso_ftest' and ( line[i]=='False' or line[i]=='FALSE' or line[i]==False):
                                replace_acceptance_criteria['check_aniso_ftest']=False

                        # search for sample criteria:
                        if header[i] in ['sample_int_n','sample_int_sigma_uT','sample_int_sigma_perc','sample_int_interval_uT','sample_int_interval_perc','sample_aniso_threshold_perc','sample_int_n_outlier_check',\
                                         'specimen_int_max_slope_diff','specimen_int_BS_68_uT','specimen_int_BS_95_uT','specimen_int_BS_68_perc','specimen_int_BS_95_perc']:
                            try:
                                replace_acceptance_criteria[header[i]]=float(line[i])
                            except:
                                pass
                        if header[i] == "sample_int_sigma":
                            try:
                                replace_acceptance_criteria['sample_int_sigma']=float(line[i])*1e-6
                                replace_acceptance_criteria['sample_int_sigma_uT']=float(line[i])*1e6
                            except:
                                pass
                        if header[i] in ["sample_int_bs_par","sample_int_bs","sample_int_stdev_opt"]:
                            if line[i]==True or line[i] in ["True","TRUE","1"]:
                                replace_acceptance_criteria[header[i]]=True
                                
            if  replace_acceptance_criteria["sample_int_bs_par"]==False and replace_acceptance_criteria["sample_int_bs"]==False and replace_acceptance_criteria["sample_int_stdev_opt"]==False:
                replace_acceptance_criteria["sample_int_stdev_opt"]=True
            
            fin.close()
            return(replace_acceptance_criteria)
        
        #except:
        #    print("-W- Cant read Criteria file from path\n")
        #    print("-I- using default criteria\n")
        except:
            return(default_acceptance_criteria)

    #----------------------------------------------------------------------

    #========================================================
    # Anistropy tensors
    #========================================================


    def calculate_anistropy_tensors(self):
        print "calling calculate_anistropy_tensors()"
        variable = """

        def tauV(T):
         
         #   gets the eigenvalues (tau) and eigenvectors (V) from matrix T
         
            print "calling tauV()"
            t,V,tr=[],[],0.
            ind1,ind2,ind3=0,1,2
            evalues,evectmps=linalg.eig(T)
            evectors=transpose(evectmps)  # to make compatible with Numeric convention
            for tau in evalues:
                tr+=tau
            if tr!=0:
                for i in range(3):
                    evalues[i]=evalues[i]/tr
            else:
                return t,V
        # sort evalues,evectors
            t1,t2,t3=0.,0.,1.
            for k in range(3):
                if evalues[k] > t1: 
                    t1,ind1=evalues[k],k 
                if evalues[k] < t3: 
                    t3,ind3=evalues[k],k 
            for k in range(3):
                if evalues[k] != t1 and evalues[k] != t3: 
                    t2,ind2=evalues[k],k
            V.append(evectors[ind1])
            V.append(evectors[ind2])
            V.append(evectors[ind3])
            t.append(t1)
            t.append(t2)
            t.append(t3)
            return t,V
        
        #def main():

            
        def calculate_aniso_parameters(B,K):
            print "calling calculate aniso_parameters()"
            aniso_parameters={}
            S_bs=dot(B,K)
            
            # normalize by trace
            trace=S_bs[0]+S_bs[1]+S_bs[2]
            S_bs=S_bs/trace
            s1,s2,s3,s4,s5,s6=S_bs[0],S_bs[1],S_bs[2],S_bs[3],S_bs[4],S_bs[5]
            s_matrix=[[s1,s4,s6],[s4,s2,s5],[s6,s5,s3]]
            
            # calculate eigen vector,
            t,evectors=eig(s_matrix)
            # sort vectors
            t=list(t)
            t1=max(t)
            ix_1=t.index(t1)
            t3=min(t)
            ix_3=t.index(t3)
            for tt in range(3):
                if t[tt]!=t1 and t[tt]!=t3:
                    t2=t[tt]
                    ix_2=t.index(t2)
                    
            v1=[evectors[0][ix_1],evectors[1][ix_1],evectors[2][ix_1]]
            v2=[evectors[0][ix_2],evectors[1][ix_2],evectors[2][ix_2]]
            v3=[evectors[0][ix_3],evectors[1][ix_3],evectors[2][ix_3]]


            DIR_v1=self.cart2dir(v1)
            DIR_v2=self.cart2dir(v2)
            DIR_v3=self.cart2dir(v3)

                               
            aniso_parameters['anisotropy_s1']="%f"%s1
            aniso_parameters['anisotropy_s2']="%f"%s2
            aniso_parameters['anisotropy_s3']="%f"%s3
            aniso_parameters['anisotropy_s4']="%f"%s4
            aniso_parameters['anisotropy_s5']="%f"%s5
            aniso_parameters['anisotropy_s6']="%f"%s6
            aniso_parameters['anisotropy_degree']="%f"%(t1/t3)
            aniso_parameters['anisotropy_t1']="%f"%t1
            aniso_parameters['anisotropy_t2']="%f"%t2
            aniso_parameters['anisotropy_t3']="%f"%t3
            aniso_parameters['anisotropy_v1_dec']="%.1f"%DIR_v1[0]
            aniso_parameters['anisotropy_v1_inc']="%.1f"%DIR_v1[1]
            aniso_parameters['anisotropy_v2_dec']="%.1f"%DIR_v2[0]
            aniso_parameters['anisotropy_v2_inc']="%.1f"%DIR_v2[1]
            aniso_parameters['anisotropy_v3_dec']="%.1f"%DIR_v3[0]
            aniso_parameters['anisotropy_v3_inc']="%.1f"%DIR_v3[1]

            # modified from pmagpy:
            if len(K)/3==9 or len(K)/3==6 or len(K)/3==15:
                n_pos=len(K)/3
                tmpH = Matrices[n_pos]['tmpH']
                a=s_matrix
                S=0.
                comp=zeros((n_pos*3),'f')
                for i in range(n_pos):
                    for j in range(3):
                        index=i*3+j
                        compare=a[j][0]*tmpH[i][0]+a[j][1]*tmpH[i][1]+a[j][2]*tmpH[i][2]
                        comp[index]=compare
                for i in range(n_pos*3):
                    d=K[i]/trace - comp[i] # del values
                    S+=d*d
                nf=float(n_pos*3-6) # number of degrees of freedom
                if S >0: 
                    sigma=math.sqrt(S/nf)
                hpars=self.dohext(nf,sigma,[s1,s2,s3,s4,s5,s6])
                
                aniso_parameters['anisotropy_sigma']="%f"%sigma
                aniso_parameters['anisotropy_ftest']="%f"%hpars["F"]
                aniso_parameters['anisotropy_ftest12']="%f"%hpars["F12"]
                aniso_parameters['anisotropy_ftest23']="%f"%hpars["F23"]
                aniso_parameters['result_description']="Critical F: %s"%(hpars['F_crit'])
                aniso_parameters['anisotropy_F_crit']="%f"%float(hpars['F_crit'])
                aniso_parameters['anisotropy_n']=n_pos
                
            return(aniso_parameters)



        
        aniso_logfile=open(self.WD+"/rmag_anisotropy.log",'w')

        aniso_logfile.write("------------------------\n")
        aniso_logfile.write( "-I- Start rmag anisrotropy script\n")
        aniso_logfile.write( "------------------------\n")



        #-----------------------------------
        # Prepare rmag_anisotropy.txt file for writing
        #-----------------------------------

        rmag_anisotropy_file =open(self.WD+"/rmag_anisotropy.txt",'w')
        rmag_anisotropy_file.write("tab\trmag_anisotropy\n")

        rmag_results_file =open(self.WD+"/rmag_results.txt",'w')
        rmag_results_file.write("tab\trmag_results\n")
        
        rmag_anistropy_header=['er_specimen_name','er_sample_name','er_site_name','anisotropy_type','anisotropy_n','anisotropy_description','anisotropy_s1','anisotropy_s2','anisotropy_s3','anisotropy_s4','anisotropy_s5','anisotropy_s6','anisotropy_sigma','anisotropy_alt','magic_experiment_names','magic_method_codes','rmag_anisotropy_name']

        String=""
        for i in range (len(rmag_anistropy_header)):
            String=String+rmag_anistropy_header[i]+'\t'
        rmag_anisotropy_file.write(String[:-1]+"\n")
        


        rmag_results_header=['er_specimen_names','er_sample_names','er_site_names','anisotropy_type','magic_method_codes','magic_experiment_names','result_description','anisotropy_t1','anisotropy_t2','anisotropy_t3','anisotropy_ftest','anisotropy_ftest12','anisotropy_ftest23',\
                             'anisotropy_v1_dec','anisotropy_v1_inc','anisotropy_v2_dec','anisotropy_v2_inc','anisotropy_v3_dec','anisotropy_v3_inc']


        String=""
        for i in range (len(rmag_results_header)):
            String=String+rmag_results_header[i]+'\t'
        rmag_results_file.write(String[:-1]+"\n")

        #-----------------------------------
        # Matrices definitions:
        # A design matrix
        # B dot(inv(dot(A.transpose(),A)),A.transpose())
        # tmpH is used for sigma calculation (9,15 measurements only)
        # 
        #  Anisotropy tensor:
        #
        # |Mx|   |s1 s4 s6|   |Bx|
        # |My| = |s4 s2 s5| . |By|
        # |Mz|   |s6 s5 s3|   |Bz|
        #
        # A matrix (measurement matrix):
        # Each mesurement yields three lines in "A" matrix
        #
        # |Mi  |   |Bx 0  0   By  0   Bz|   |s1|
        # |Mi+1| = |0  By 0   Bx  Bz  0 | . |s2|
        # |Mi+2|   |0  0  Bz  0   By  Bx|   |s3|
        #                                   |s4|
        #                                   |s5|
        #
        #-----------------------------------

        Matrices={}
        
        for n_pos in [6,9,15]:

            Matrices[n_pos]={}
            
            A=zeros((n_pos*3,6),'f')

            if n_pos==6:
                positions=[[0.,0.,1.],[90.,0.,1.],[0.,90.,1.],\
                     [180.,0.,1.],[270.,0.,1.],[0.,-90.,1.]]

            if n_pos==15:
                positions=[[315.,0.,1.],[225.,0.,1.],[180.,0.,1.],[135.,0.,1.],[45.,0.,1.],\
                     [90.,-45.,1.],[270.,-45.,1.],[270.,0.,1.],[270.,45.,1.],[90.,45.,1.],\
                     [180.,45.,1.],[180.,-45.,1.],[0.,-90.,1.],[0,-45.,1.],[0,45.,1.]]
            if n_pos==9:
                positions=[[315.,0.,1.],[225.,0.,1.],[180.,0.,1.],\
                     [90.,-45.,1.],[270.,-45.,1.],[270.,0.,1.],\
                     [180.,45.,1.],[180.,-45.,1.],[0.,-90.,1.]]

            
            tmpH=zeros((n_pos,3),'f') # define tmpH
            for i in range(len(positions)):
                CART=self.dir2cart(positions[i])
                a=CART[0];b=CART[1];c=CART[2]
                A[3*i][0]=a
                A[3*i][3]=b
                A[3*i][5]=c

                A[3*i+1][1]=b
                A[3*i+1][3]=a
                A[3*i+1][4]=c

                A[3*i+2][2]=c
                A[3*i+2][4]=b
                A[3*i+2][5]=a
                
                tmpH[i][0]=CART[0]
                tmpH[i][1]=CART[1]
                tmpH[i][2]=CART[2]

            B=dot(inv(dot(A.transpose(),A)),A.transpose())

            Matrices[n_pos]['A']=A
            Matrices[n_pos]['B']=B
            Matrices[n_pos]['tmpH']=tmpH

        #==================================================================================

        Data_anisotropy={}                
        specimens=self.Data.keys()
        specimens.sort()
        for specimen in specimens:

            if 'atrmblock' in self.Data[specimen].keys():
                
                #-----------------------------------
                # aTRM 6 positions
                #-----------------------------------
                    
                #print "-I- Start calculating ATRM tensor for "
                atrmblock=self.Data[specimen]['atrmblock']
                trmblock=self.Data[specimen]['trmblock']
                zijdblock=self.Data[specimen]['zijdblock']
                if len(atrmblock)<6:
                    aniso_logfile.write("-W- specimen %s has not enough measurementf for ATRM calculation\n"%specimen)
                    continue
                
                B=Matrices[6]['B']
                                    
                Reject_specimen = False

                # The zero field step is a "baseline"
                # and the atrm measurements are substructed from the baseline
                # if there is a zero field is in the atrm block: then use this measurement as a baseline
                # if not, the program searches for the zero-field step in the zijderveld block. 
                # the baseline is the average of all the zero field steps in the same temperature (in case there is more than one)

                # Search the baseline in the ATRM measurement
                # Search the alteration check in the ATRM measurement
                # If there is more than one baseline measurements then avrage all measurements
                
                baseline=""
                Alteration_check=""
                Alteration_check_index=""
                baselines=[]

                # search for baseline in atrm blocks
                for rec in atrmblock:
                    dec=float(rec['measurement_dec'])
                    inc=float(rec['measurement_inc'])
                    moment=float(rec['measurement_magn_moment'])
                    # find the temperature of the atrm
                    if float(rec['treatment_dc_field'])!=0 and float(rec['treatment_temp'])!=273:
                        atrm_temperature=float(rec['treatment_temp'])
                    # find baseline
                    if float(rec['treatment_dc_field'])==0 and float(rec['treatment_temp'])!=273:
                        baselines.append(array(self.dir2cart([dec,inc,moment])))
                    # Find alteration check
                    #print rec['measurement_number']
                
                if len(baselines)!=0:
                    aniso_logfile.write( "-I- found ATRM baseline for specimen %s\n"%specimen)
                    
                else:
                    if len(zijdblock)!=0 :
                        for rec in zijdblock:
                            zij_temp=rec[0]
                            #print rec
                            if zij_temp==atrm_temperature:
                                dec=float(rec[1])
                                inc=float(rec[2])
                                moment=float(rec[3])
                                baselines.append(array(self.dir2cart([dec,inc,moment])))
                                aniso_logfile.write( "-I- Found %i ATRM baselines for specimen %s in Zijderveld block. Averaging measurements\n"%(len(baselines),specimen))
                if  len(baselines)==0:
                    baseline=zeros(3,'f')
                    aniso_logfile.write( "-I- No aTRM baseline for specimen %s\n"%specimen)
                else:
                    baselines=array(baselines)
                    baseline=array([mean(baselines[:,0]),mean(baselines[:,1]),mean(baselines[:,2])])                                 
                           
                # sort measurements
                
                M=zeros([6,3],'f')
                
                for rec in atrmblock:

                    dec=float(rec['measurement_dec'])
                    inc=float(rec['measurement_inc'])
                    moment=float(rec['measurement_magn_moment'])
                    CART=array(self.dir2cart([dec,inc,moment]))-baseline
                    
                    if float(rec['treatment_dc_field'])==0: # Ignore zero field steps
                        continue
                    if  "LT-PTRM-I" in rec['magic_method_codes'].split(":"): #  alteration check
                        Alteration_check=CART
                        Alteration_check_dc_field_phi=float(rec['treatment_dc_field_phi'])
                        Alteration_check_dc_field_theta=float(rec['treatment_dc_field_theta'])
                        if Alteration_check_dc_field_phi==0 and Alteration_check_dc_field_theta==0 :
                            Alteration_check_index=0
                        if Alteration_check_dc_field_phi==90 and Alteration_check_dc_field_theta==0 :
                            Alteration_check_index=1
                        if Alteration_check_dc_field_phi==0 and Alteration_check_dc_field_theta==90 :
                            Alteration_check_index=2
                        if Alteration_check_dc_field_phi==180 and Alteration_check_dc_field_theta==0 :
                            Alteration_check_index=3
                        if Alteration_check_dc_field_phi==270 and Alteration_check_dc_field_theta==0 :
                            Alteration_check_index=4
                        if Alteration_check_dc_field_phi==0 and Alteration_check_dc_field_theta==-90 :
                            Alteration_check_index=5
                        aniso_logfile.write(  "-I- found alteration check  for specimen %s\n"%specimen)
                        continue
                    
                    treatment_dc_field_phi=float(rec['treatment_dc_field_phi'])
                    treatment_dc_field_theta=float(rec['treatment_dc_field_theta'])
                    treatment_dc_field=float(rec['treatment_dc_field'])
                    
                    #+x, M[0]
                    if treatment_dc_field_phi==0 and treatment_dc_field_theta==0 :
                        M[0]=CART
                    #+Y , M[1]
                    if treatment_dc_field_phi==90 and treatment_dc_field_theta==0 :
                        M[1]=CART
                    #+Z , M[2]
                    if treatment_dc_field_phi==0 and treatment_dc_field_theta==90 :
                        M[2]=CART
                    #-x, M[3]
                    if treatment_dc_field_phi==180 and treatment_dc_field_theta==0 :
                        M[3]=CART
                    #-Y , M[4]
                    if treatment_dc_field_phi==270 and treatment_dc_field_theta==0 :
                        M[4]=CART
                    #-Z , M[5]
                    if treatment_dc_field_phi==0 and treatment_dc_field_theta==-90 :
                        M[5]=CART
            
                # check if at least one measurement in missing
                for i in range(len(M)):
                    if M[i][0]==0 and M[i][1]==0 and M[i][2]==0: 
                        aniso_logfile.write( "-E- ERROR: missing atrm data for specimen %s\n"%(specimen))
                        Reject_specimen=True

                # alteration check        

                anisotropy_alt=0
                if Alteration_check!="":
                    for i in range(len(M)):
                        if Alteration_check_index==i:
                            M_1=sqrt(sum((array(M[i])**2)))
                            M_2=sqrt(sum(Alteration_check**2))
                            diff=abs(M_1-M_2)
                            diff_ratio=diff/max(M_1,M_2)
                            diff_ratio_perc=100*diff_ratio
                            if diff_ratio_perc > anisotropy_alt:
                                anisotropy_alt=diff_ratio_perc
                else:
                    aniso_logfile.write( "-W- Warning: no alteration check for specimen %s \n "%specimen )

                # Check for maximum difference in anti parallel directions.
                # if the difference between the two measurements is more than maximum_diff
                # The specimen is rejected
                
                # i.e. +x versus -x, +y versus -y, etc.s

                for i in range(3):
                    M_1=sqrt(sum(array(M[i])**2))
                    M_2=sqrt(sum(array(M[i+3])**2))
                    
                    diff=abs(M_1-M_2)
                    diff_ratio=diff/max(M_1,M_2)
                    diff_ratio_perc=100*diff_ratio
                    
                    if diff_ratio_perc>anisotropy_alt:
                        anisotropy_alt=diff_ratio_perc
                        
                if not Reject_specimen:
                
                    # K vector (18 elements, M1[x], M1[y], M1[z], ... etc.) 
                    K=zeros(18,'f')
                    K[0],K[1],K[2]=M[0][0],M[0][1],M[0][2]
                    K[3],K[4],K[5]=M[1][0],M[1][1],M[1][2]
                    K[6],K[7],K[8]=M[2][0],M[2][1],M[2][2]
                    K[9],K[10],K[11]=M[3][0],M[3][1],M[3][2]
                    K[12],K[13],K[14]=M[4][0],M[4][1],M[4][2]
                    K[15],K[16],K[17]=M[5][0],M[5][1],M[5][2]

                    if specimen not in Data_anisotropy.keys():
                        Data_anisotropy[specimen]={}
                    aniso_parameters=calculate_aniso_parameters(B,K)
                    Data_anisotropy[specimen]['ATRM']=aniso_parameters
                    Data_anisotropy[specimen]['ATRM']['anisotropy_alt']="%.2f"%anisotropy_alt               
                    Data_anisotropy[specimen]['ATRM']['anisotropy_type']="ATRM"
                    Data_anisotropy[specimen]['ATRM']['er_sample_name']=atrmblock[0]['er_sample_name']
                    Data_anisotropy[specimen]['ATRM']['er_specimen_name']=specimen
                    Data_anisotropy[specimen]['ATRM']['er_site_name']=atrmblock[0]['er_site_name']
                    Data_anisotropy[specimen]['ATRM']['anisotropy_description']='Hext statistics adapted to ATRM'
                    Data_anisotropy[specimen]['ATRM']['magic_experiment_names']=specimen+";ATRM"
                    Data_anisotropy[specimen]['ATRM']['magic_method_codes']="LP-AN-TRM:AE-H"
                    Data_anisotropy[specimen]['ATRM']['rmag_anisotropy_name']=specimen


            if 'aarmblock' in self.Data[specimen].keys():    

                #-----------------------------------
                # AARM - 6, 9 or 15 positions
                #-----------------------------------
                    
                aniso_logfile.write( "-I- Start calculating AARM tensors specimen %s\n"%specimen)

                aarmblock=self.Data[specimen]['aarmblock']
                if len(aarmblock)<12:
                    aniso_logfile.write( "-W- WARNING: not enough aarm measurement for specimen %s\n"%specimen)
                    continue
                elif len(aarmblock)==12:
                    n_pos=6
                    B=Matrices[6]['B']
                    M=zeros([6,3],'f')
                elif len(aarmblock)==18:
                    n_pos=9
                    B=Matrices[9]['B']
                    M=zeros([9,3],'f')
                # 15 positions
                elif len(aarmblock)==30:
                    n_pos=15
                    B=Matrices[15]['B']
                    M=zeros([15,3],'f')
                else:
                    aniso_logfile.write( "-E- ERROR: number of measurements in aarm block is incorrect sample %s\n"%specimen)
                    continue
                    
                Reject_specimen = False

                for i in range(n_pos):
                    for rec in aarmblock:
                        if float(rec['measurement_number'])==i*2+1:
                            dec=float(rec['measurement_dec'])
                            inc=float(rec['measurement_inc'])
                            moment=float(rec['measurement_magn_moment'])                    
                            M_baseline=array(self.dir2cart([dec,inc,moment]))
                            
                        if float(rec['measurement_number'])==i*2+2:
                            dec=float(rec['measurement_dec'])
                            inc=float(rec['measurement_inc'])
                            moment=float(rec['measurement_magn_moment'])                    
                            M_arm=array(self.dir2cart([dec,inc,moment]))
                    M[i]=M_arm-M_baseline

                    
                K=zeros(3*n_pos,'f')
                for i in range(n_pos):
                    K[i*3]=M[i][0]
                    K[i*3+1]=M[i][1]
                    K[i*3+2]=M[i][2]            

                if specimen not in Data_anisotropy.keys():
                    Data_anisotropy[specimen]={}
                aniso_parameters=calculate_aniso_parameters(B,K)
                Data_anisotropy[specimen]['AARM']=aniso_parameters
                Data_anisotropy[specimen]['AARM']['anisotropy_alt']=""               
                Data_anisotropy[specimen]['AARM']['anisotropy_type']="AARM"
                Data_anisotropy[specimen]['AARM']['er_sample_name']=aarmblock[0]['er_sample_name']
                Data_anisotropy[specimen]['AARM']['er_site_name']=aarmblock[0]['er_site_name']
                Data_anisotropy[specimen]['AARM']['er_specimen_name']=specimen
                Data_anisotropy[specimen]['AARM']['anisotropy_description']='Hext statistics adapted to AARM'
                Data_anisotropy[specimen]['AARM']['magic_experiment_names']=specimen+";AARM"
                Data_anisotropy[specimen]['AARM']['magic_method_codes']="LP-AN-ARM:AE-H"
                Data_anisotropy[specimen]['AARM']['rmag_anisotropy_name']=specimen
                
    
        #-----------------------------------   

        specimens=Data_anisotropy.keys()
        specimens.sort

        # remove previous anistropy data, and replace with the new one:
        s_list=self.Data.keys()
        for sp in s_list:
            if 'AniSpec' in self.Data[sp].keys():
                del  self.Data[sp]['AniSpec']
        for specimen in specimens:
            # if both AARM and ATRM axist prefer the AARM !!
            if 'AARM' in Data_anisotropy[specimen].keys():
                TYPES=['AARM']
            if 'ATRM' in Data_anisotropy[specimen].keys():
                TYPES=['ATRM']
            if  'AARM' in Data_anisotropy[specimen].keys() and 'ATRM' in Data_anisotropy[specimen].keys():
                TYPES=['ATRM','AARM']
                aniso_logfile.write( "-W- WARNING: both aarm and atrm data exist for specimen %s. using AARM by default. If you prefer using one of them, delete the other!\n"%specimen)
            for TYPE in TYPES:
                String=""
                for i in range (len(rmag_anistropy_header)):
                    try:
                        String=String+Data_anisotropy[specimen][TYPE][rmag_anistropy_header[i]]+'\t'
                    except:
                        String=String+"%f"%(Data_anisotropy[specimen][TYPE][rmag_anistropy_header[i]])+'\t'
                rmag_anisotropy_file.write(String[:-1]+"\n")

                String=""
                Data_anisotropy[specimen][TYPE]['er_specimen_names']=Data_anisotropy[specimen][TYPE]['er_specimen_name']
                Data_anisotropy[specimen][TYPE]['er_sample_names']=Data_anisotropy[specimen][TYPE]['er_sample_name']
                Data_anisotropy[specimen][TYPE]['er_site_names']=Data_anisotropy[specimen][TYPE]['er_site_name']
                for i in range (len(rmag_results_header)):
                    try:
                        String=String+Data_anisotropy[specimen][TYPE][rmag_results_header[i]]+'\t'
                    except:
                        String=String+"%f"%(Data_anisotropy[specimen][TYPE][rmag_results_header[i]])+'\t'
                rmag_results_file.write(String[:-1]+"\n")

                if 'AniSpec' not in self.Data[specimen]:
                    self.Data[specimen]['AniSpec']={}
                self.Data[specimen]['AniSpec'][TYPE]=Data_anisotropy[specimen][TYPE]
        rmag_anisotropy_file.close()

    #==================================================        

    

        def find_sample_min_max_interpretation (Intensities,acceptance_criteria):
          print "calling find_sample_min_max_interpretation()"
          # Find the minimum and maximum acceptable sample mean.

          # make a new dictionary named "tmp_Intensities" with all grade A interpretation sorted. 
          tmp_Intensities={}
          Acceptable_sample_min,Acceptable_sample_max="",""
          for this_specimen in Intensities.keys():
            B_list=[B  for B in Intensities[this_specimen]]
            if len(B_list)>0:
                B_list.sort()
                tmp_Intensities[this_specimen]=B_list

          # find the minmum acceptable values
          while len(tmp_Intensities.keys())>=float(acceptance_criteria["sample_int_n"]):
              B_tmp=[]
              B_tmp_min=1e10
              for specimen in tmp_Intensities.keys():
                  B_tmp.append(min(tmp_Intensities[specimen]))
                  if min(tmp_Intensities[specimen])<B_tmp_min:
                      specimen_to_remove=specimen
                      B_tmp_min=min(tmp_Intensities[specimen])
              if std(B_tmp,ddof=1)<=acceptance_criteria["sample_int_sigma_uT"] or 100*(std(B_tmp,ddof=1)/mean(B_tmp))<=acceptance_criteria["sample_int_sigma_perc"]:
                  Acceptable_sample_min=mean(B_tmp)
                  #print "min value,std,",mean(B_tmp),std(B_tmp),100*(std(B_tmp)/mean(B_tmp))
                  break
              else:
                  tmp_Intensities[specimen_to_remove].remove(B_tmp_min)
                  if len(tmp_Intensities[specimen_to_remove])==0:
                      break
                      
          tmp_Intensities={}
          for this_specimen in Intensities.keys():
            B_list=[B  for B in Intensities[this_specimen]]
            if len(B_list)>0:
                B_list.sort()
                tmp_Intensities[this_specimen]=B_list

          while len(tmp_Intensities.keys())>=float(acceptance_criteria["sample_int_n"]):
              B_tmp=[]
              B_tmp_max=0
              for specimen in tmp_Intensities.keys():
                  B_tmp.append(max(tmp_Intensities[specimen]))
                  if max(tmp_Intensities[specimen])>B_tmp_max:
                      specimen_to_remove=specimen
                      B_tmp_max=max(tmp_Intensities[specimen])
              if std(B_tmp,ddof=1)<=acceptance_criteria["sample_int_sigma_uT"] or 100*(std(B_tmp,ddof=1)/mean(B_tmp))<=acceptance_criteria["sample_int_sigma_perc"]:
                  Acceptable_sample_max=mean(B_tmp)
                  #print "max value,std,",mean(B_tmp),std(B_tmp),100*(std(B_tmp)/mean(B_tmp))

                  break
              else:
                  tmp_Intensities[specimen_to_remove].remove(B_tmp_max)
                  if len(tmp_Intensities[specimen_to_remove])<1:
                      break

          if Acceptable_sample_min=="" or Acceptable_sample_max=="":
              return(0.,0.)
          return(Acceptable_sample_min,Acceptable_sample_max) 

        ############
        # End function definitions
        ############

# LJ  PROBABLY DO NOT NEED WHAT IS BELOW....
    
        start_time=time.time()
        #------------------------------------------------
        # Clean work directory
        #------------------------------------------------

        self.write_acceptance_criteria_to_file()
        try:
            shutil.rmtree(self.WD+"/thellier_interpreter")
        except:
            pass
    
        try:
            os.mkdir(self.WD+"/thellier_interpreter")
        except:
            pass

        parameters_with_upper_bounds= ['specimen_gmax','specimen_b_beta','specimen_dang','specimen_drats','specimen_int_mad','specimen_md']
        parameters_with_lower_bounds= ['specimen_int_n','specimen_int_ptrm_n','specimen_f','specimen_fvds','specimen_frac']
        accept_specimen_keys=['specimen_int_n','specimen_int_ptrm_n','specimen_f','specimen_fvds','specimen_frac','specimen_gmax','specimen_b_beta','specimen_dang','specimen_drats','specimen_int_mad','specimen_md','specimen_g','specimen_q']

        #------------------------------------------------
        # Intialize interpreter output files:
        # Prepare header for "Thellier_auto_interpretation.all.txt" 
        # All the acceptable interpretation are saved in this file
        #------------------------------------------------

        # log file
        thellier_interpreter_log=open(self.WD+"/"+"/thellier_interpreter//thellier_interpreter.log",'w')
        thellier_interpreter_log.write("-I- Start auto interpreter\n")

        # "all grade A interpretation
        thellier_interpreter_all=open(self.WD+"/thellier_interpreter/thellier_interpreter_all.txt",'w')
        thellier_interpreter_all.write("tab\tpmag_specimens\n")
        String="er_specimen_name\tmeasurement_step_min\tmeasurement_step_max\tspecimen_lab_field_dc_uT\tspecimen_int_corr_anisotropy\tspecimen_int_corr_nlt\tspecimen_int_corr_cooling_rate\tspecimen_int_uT\t"
        for key in accept_specimen_keys:
            String=String+key+"\t"
        String=String[:-1]+"\n"
        thellier_interpreter_all.write(String)

        #specimen_bound
        Fout_specimens_bounds=open(self.WD+"/thellier_interpreter/thellier_interpreter_specimens_bounds.txt",'w')
        String="Selection criteria:\n"
        for key in accept_specimen_keys:
                String=String+key+"\t"
        Fout_specimens_bounds.write(String[:-1]+"\n")
        String=""
        for key in accept_specimen_keys:
            if key!= "specimen_frac":
                String=String+"%.2f\t"%self.accept_new_parameters[key]
            else:
                String=String+"%s\t"%self.accept_new_parameters[key]                
        Fout_specimens_bounds.write(String[:-1]+"\n")
        
        Fout_specimens_bounds.write("--------------------------------\n")
        Fout_specimens_bounds.write("er_sample_name\ter_specimen_name\tspecimen_int_corr_anisotropy\tAnisotropy_code\tspecimen_int_corr_nlt\tspecimen_int_corr_cooling_rate\tspecimen_lab_field_dc_uT\tspecimen_int_min_uT\tspecimen_int_max_uT\tWARNING\n")


        criteria_string="Selection criteria:\n"
        for key in self.accept_new_parameters.keys():
            if "sample" in key:
                criteria_string=criteria_string+key+"\t"
        for key in accept_specimen_keys:
            if "specimen" in key:
                criteria_string=criteria_string+key+"\t"
        criteria_string=criteria_string[:-1]+"\n"
        for key in self.accept_new_parameters.keys():
            if "sample" in key:
                try:
                    criteria_string=criteria_string+"%.2f"%self.accept_new_parameters[key]+"\t"
                except:
                    criteria_string=criteria_string+"%s"%self.accept_new_parameters[key]+"\t"                
        for key in accept_specimen_keys:
            if "specimen" in key:
                try:
                    criteria_string=criteria_string+"%.2f"%self.accept_new_parameters[key]+"\t"
                except:
                    criteria_string=criteria_string+"%s"%self.accept_new_parameters[key]+"\t"
        criteria_string=criteria_string[:-1]+"\n"
        criteria_string=criteria_string+"---------------------------------\n"

        # STDEV-OPT output files
        if self.accept_new_parameters['sample_int_stdev_opt']:
            Fout_STDEV_OPT_redo=open(self.WD+"/thellier_interpreter/thellier_interpreter_STDEV-OPT_redo",'w')

            Fout_STDEV_OPT_specimens=open(self.WD+"/thellier_interpreter/thellier_interpreter_STDEV-OPT_specimens.txt",'w')
            Fout_STDEV_OPT_specimens.write("tab\tpmag_specimens\n")
            String="er_sample_name\ter_specimen_name\tspecimen_int_uT\tmeasurement_step_min\tmeasurement_step_min\tspecimen_lab_field_dc\tAnisotropy_correction_factor\tNLT_correction_factor\tCooling_rate_correction_factor\t"
            for key in accept_specimen_keys:
                String=String+key+"\t"        
            Fout_STDEV_OPT_specimens.write(String[:-1]+"\n")

            Fout_STDEV_OPT_samples=open(self.WD+"/thellier_interpreter/thellier_interpreter_STDEV-OPT_samples.txt",'w')
            Fout_STDEV_OPT_samples.write(criteria_string)
            Fout_STDEV_OPT_samples.write("er_sample_name\tsample_int_n\tsample_int_uT\tsample_int_sigma_uT\tsample_int_sigma_perc\tsample_int_interval_uT\tsample_int_interval_perc\tWarning\n")
        # simple bootstrap output files
        

                
 
        if self.accept_new_parameters['sample_int_bs']:
           Fout_BS_samples=open(self.WD+"/thellier_interpreter/thellier_interpreter_BS_samples.txt",'w')
           Fout_BS_samples.write(criteria_string)
           #Fout_BS_samples.write("---------------------------------\n")
           Fout_BS_samples.write("er_sample_name\tsample_int_n\tsample_int_uT\tsample_int_68_low\tsample_int_68_high\tsample_int_95_low\tsample_int_95_high\tsample_int_sigma_uT\tsample_int_sigma_perc\tWARNING\n")
        # parameteric bootstrap output files

        if self.accept_new_parameters['sample_int_bs_par']:
           Fout_BS_PAR_samples=open(self.WD+"/thellier_interpreter/thellier_interpreter_BS-PAR_samples.txt",'w')
           Fout_BS_PAR_samples.write(criteria_string) 
           #Fout_BS_PAR_samples.write("---------------------------------\n")
           Fout_BS_PAR_samples.write("er_sample_name\tsample_int_n\tsample_int_uT\tsample_int_68_low\tsample_int_68_high\tsample_int_95_low\tsample_int_95_high\tsample_int_sigma_uT\tsample_int_sigma_perc\tWARNING\n")
           
        thellier_interpreter_log.write("-I- using paleointenisty statistics:\n")
        for key in [key for key in self.accept_new_parameters.keys() if "sample" in key]:
            try:
                thellier_interpreter_log.write("-I- %s=%.2f\n"%(key,self.accept_new_parameters[key]))
            except:
                thellier_interpreter_log.write("-I- %s=%s\n"%(key,self.accept_new_parameters[key]))
                                            
        for key in [key for key in self.accept_new_parameters.keys() if "specimen" in key]:
            try:
                thellier_interpreter_log.write("-I- %s=%.2f\n"%(key,self.accept_new_parameters[key]))
            except:
                thellier_interpreter_log.write("-I- %s=%s\n"%(key,self.accept_new_parameters[key]))
               
                                  

    
        #------------------------------------------------
    
        busy_frame=wx.BusyInfo("Running Thellier auto interpreter\n It may take several minutes depending on the number of specimens ...", self)

        specimens_list=self.Data.keys()
        specimens_list.sort()
        thellier_interpreter_log.write("-I- Found %i specimens\n"%(len(specimens_list)))

        #try:
        All_grade_A_Recs={}
        for s in specimens_list:
            thellier_interpreter_log.write("-I- doing now specimen %s\n"%s)
            self.Data[s]['pars']={}
            self.Data[s]['pars']['lab_dc_field']=self.Data[s]['lab_dc_field']
            self.Data[s]['pars']['er_specimen_name']=s
            self.Data[s]['pars']['er_sample_name']=self.Data_hierarchy['specimens'][s]
            temperatures=self.Data[s]['t_Arai']
            
            # check that all temperatures are in right order:
            ignore_specimen=False
            for t in range(len(temperatures)-1):
                if float(temperatures[t+1])<float(temperatures[t]):
                    thellier_interpreter_log.write("-W- Found problem in the temperature order of specimen %s. skipping specimen\n"%(s))
                    ignore_specimen=True
            if ignore_specimen:
                continue
            specimen_int_n=int(self.accept_new_parameters['specimen_int_n'])
            for tmin_i in range(len(temperatures)-specimen_int_n+1):
                for tmax_i in range(tmin_i+specimen_int_n-1,len(temperatures)):
                    tmin=temperatures[tmin_i]
                    tmax=temperatures[tmax_i]
                    pars=self.get_PI_parameters(s,tmin,tmax)

                    #-------------------------------------------------            
                    # check if pass the criteria
                    #-------------------------------------------------

                    # distinguish between upper threshold value and lower threshold value

                    if 'specimen_fail_criteria' not in pars:
                        continue
                    if len(pars['specimen_fail_criteria'])>0:
                        # Fail:
                        message_string= "-I- specimen %s (%.0f-%.0f) FAIL on: "%(s,float(pars["measurement_step_min"])-273, float(pars["measurement_step_max"])-273)
                        for parameter in pars['specimen_fail_criteria']:
                            if "scat" not in parameter:
                                message_string=message_string+parameter + "= %f,  "%pars[parameter]
                            else:
                                message_string=message_string+parameter + "= %s,  "%str(pars[parameter])
                                
                        thellier_interpreter_log.write(message_string+"\n")        
                    else:

                        # PASS:
                        message_string = "-I- specimen %s (%.0f-%.0f) PASS"%(s,float(pars["measurement_step_min"])-273, float(pars["measurement_step_max"])-273)
                        thellier_interpreter_log.write(message_string+"\n")
                        
                        #--------------------------------------------------------------
                        # Save all the grade A interpretation in thellier_interpreter_all.txt
                        #--------------------------------------------------------------

                        String=s+"\t"
                        String=String+"%.0f\t"%(float(pars["measurement_step_min"])-273.)
                        String=String+"%.0f\t"%(float(pars["measurement_step_max"])-273.)
                        String=String+"%.0f\t"%(float(pars["lab_dc_field"])*1e6)
                        if "AC_specimen_correction_factor" in pars.keys():
                           String=String+"%.2f\t"%float(pars["AC_specimen_correction_factor"])
                        else:
                           String=String+"-\t"
                        if  float(pars["NLT_specimen_correction_factor"])!=-999:
                           String=String+"%.2f\t"%float(pars["NLT_specimen_correction_factor"])
                        else:
                           String=String+"-\t"
                        if  float(pars["specimen_int_corr_cooling_rate"])!=-999 and float(pars["specimen_int_corr_cooling_rate"])!=-1 :
                           String=String+"%.2f\t"%float(pars["specimen_int_corr_cooling_rate"])
                        else:
                           String=String+"-\t"
                        Bancient=float(pars['specimen_int_uT'])
                        String=String+"%.1f\t"%(Bancient)

                        for key in accept_specimen_keys:
                           String=String+"%.2f"%(pars[key])+"\t"
                        String=String[:-1]+"\n"

                        thellier_interpreter_all.write(String)


                        #-------------------------------------------------                    
                        # save 'acceptable' (grade A) specimen interpretaion
                        #-------------------------------------------------
                        
                        if s not in All_grade_A_Recs.keys():
                           All_grade_A_Recs[s]={}
                        new_pars={}
                        for k in pars.keys():
                            new_pars[k]=pars[k]
                        TEMP="%.0f,%.0f"%(float(pars["measurement_step_min"])-273,float(pars["measurement_step_max"])-273)
                        All_grade_A_Recs[s][TEMP]=new_pars


        specimens_list=All_grade_A_Recs.keys()
        specimens_list.sort()
        Grade_A_samples={}
        Redo_data_specimens={}

        
        #--------------------------------------------------------------
        # specimens bound file
        #--------------------------------------------------------------



        for s in specimens_list:

            sample=self.Data_hierarchy['specimens'][s]
            B_lab=float(self.Data[s]['lab_dc_field'])*1e6
            B_min,B_max=1e10,0.
            NLT_factor_min,NLT_factor_max=1e10,0.
            all_B_tmp_array=[]

            for TEMP in All_grade_A_Recs[s].keys():
                pars=All_grade_A_Recs[s][TEMP]
                if "AC_anisotropy_type" in pars.keys():
                    AC_correction_factor=pars["Anisotropy_correction_factor"]
                    AC_correction_type=pars["AC_anisotropy_type"]
                    WARNING=""
                    if "AC_WARNING" in pars.keys():
                        WARNING=WARNING+pars["AC_WARNING"]
                else:
                    AC_correction_factor=1.
                    AC_correction_type="-"
                    WARNING="WARNING: No anisotropy correction"
                
                B_anc=pars['specimen_int_uT']
                    
                if B_anc< B_min:
                    B_min=B_anc
                if B_anc > B_max:
                    B_max=B_anc
                if pars["NLT_specimen_correction_factor"]!=-1:
                    NLT_f=pars['NLT_specimen_correction_factor']
                    if NLT_f< NLT_factor_min:
                        NLT_factor_min=NLT_f
                    if NLT_f > NLT_factor_max:
                        NLT_factor_max=NLT_f                

                # sort by samples
                #--------------------------------------------------------------
                
                if sample not in Grade_A_samples.keys():
                    Grade_A_samples[sample]={}
                if s not in Grade_A_samples[sample].keys() and len(All_grade_A_Recs[s])>0:
                    Grade_A_samples[sample][s]=[]
                if s not in Redo_data_specimens.keys():
                    Redo_data_specimens[s]={}

                Grade_A_samples[sample][s].append(B_anc)                
                #Redo_data_specimens[s][B_anc]=pars

            # write to specimen_bounds
            #--------------------------------------------------------------

            if pars["NLT_specimen_correction_factor"] != -1:
                NLT_factor="%.2f"%(NLT_factor_max)
            else:
                NLT_factor="-"

            if pars["specimen_int_corr_cooling_rate"] != -1 and pars["specimen_int_corr_cooling_rate"] != -999:
                CR_factor="%.2f"%(float(pars["specimen_int_corr_cooling_rate"]))
            else:
                CR_factor="-"
            if 'cooling_rate_data' in  self.Data[s].keys():
                if 'CR_correction_factor_flag' in  self.Data[s]['cooling_rate_data'].keys():
                    if self.Data[s]['cooling_rate_data']['CR_correction_factor_flag'] != "calculated":
                        if "inferred" in self.Data[s]['cooling_rate_data']['CR_correction_factor_flag']:
                            WARNING=WARNING+";"+"cooling rate correction inferred from sister specimens"
                        if "alteration" in self.Data[s]['cooling_rate_data']['CR_correction_factor_flag']:
                            WARNING=WARNING+";"+"cooling rate experiment failed alteration"
                        if "bad" in self.Data[s]['cooling_rate_data']['CR_correction_factor_flag']:
                            WARNING=WARNING+";"+"cooling rate experiment failed"
                
            if AC_correction_type =="-":
                AC_correction_factor_to_print="-"
            else:
                AC_correction_factor_to_print="%.2f"%AC_correction_factor
            
            String="%s\t%s\t%s\t%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%s\n"\
                    %(sample,s,AC_correction_factor_to_print,AC_correction_type,NLT_factor,CR_factor,B_lab,B_min,B_max,WARNING)
            Fout_specimens_bounds.write(String)


        #--------------------------------------------------------------
        # Find the STDEV-OPT 'best mean':
        # the interprettaions that give
        # the minimum standrad deviation.
        #
        #--------------------------------------------------------------

        # Sort all grade A interpretation

        samples=Grade_A_samples.keys()
        samples.sort()

        # show all interpretation after sone with the interpreter. 
        # first delete all previous interpretation
        for sp in self.Data.keys():
            del self.Data[sp]['pars']
            self.Data[sp]['pars']={}
            self.Data[sp]['pars']['lab_dc_field']=self.Data[sp]['lab_dc_field']
            self.Data[sp]['pars']['er_specimen_name']=self.Data[sp]['er_specimen_name']   
            self.Data[sp]['pars']['er_sample_name']=self.Data[sp]['er_sample_name']
        self.Data_samples={}
        interpreter_redo={}

        

        for sample in samples:

            #--------------------------------------------------------------
            # check for anistropy issue:
            # If the average anisotropy correction in the sample is > 10%,
            # and there are enough good specimens with  anisotropy correction to pass sample's criteria
            # then dont use the uncorrected specimens for sample's calculation. 
            #--------------------------------------------------------------

            if self.accept_new_parameters['sample_aniso_threshold_perc'] != "" and float(self.accept_new_parameters['sample_aniso_threshold_perc'])<100:
                if len(Grade_A_samples[sample].keys())>self.accept_new_parameters['sample_int_n']:
                    aniso_corrections=[]
                    for specimen in Grade_A_samples[sample].keys():
                        AC_correction_factor=0
                        for k in All_grade_A_Recs[specimen].keys():
                            pars=All_grade_A_Recs[specimen][k]
                            if "AC_anisotropy_type" in pars.keys():
                                if "AC_WARNING" in pars.keys():
                                    if "TRM" in pars["AC_WARNING"] and  pars["AC_anisotropy_type"]== "ATRM" \
                                       or "ARM" in pars["AC_WARNING"] and  pars["AC_anisotropy_type"]== "AARM":
                                        continue
                                    AC_correction_factor=max(AC_correction_factor,pars["Anisotropy_correction_factor"])
                        if AC_correction_factor!=0:
                            aniso_corrections.append(abs(1.-float(AC_correction_factor)))
                    if aniso_corrections!=[]:
                        thellier_interpreter_log.write("sample %s have anisotropy factor mean of %f\n"%(sample,mean(aniso_corrections)))

                    if mean(aniso_corrections) > float(self.accept_new_parameters['sample_aniso_threshold_perc'])/100 : # 0.10:
                        tmp_Grade_A_samples=copy.deepcopy(Grade_A_samples)
                        warning_messeage=""
                        WARNING_tmp=""
                        #print "sample %s have anisotropy factor mean of %f"%(sample,mean(aniso_corrections))
                        for specimen in Grade_A_samples[sample].keys():
                            ignore_specimen=False
                            intenstities=All_grade_A_Recs[specimen].keys()
                            pars=All_grade_A_Recs[specimen][intenstities[0]]
                            if "AC_anisotropy_type" not in pars.keys():
                                ignore_specimen=True
                            elif "AC_WARNING" in pars.keys():
                                if "alteration check" in pars["AC_WARNING"]:
                                    if "TRM" in pars["AC_WARNING"] and  pars["AC_anisotropy_type"]== "ATRM" \
                                       or "ARM" in pars["AC_WARNING"] and  pars["AC_anisotropy_type"]== "AARM":
                                        ignore_specimen=True
                            if ignore_specimen: 
                                warning_messeage = warning_messeage + "-W- WARNING: specimen %s is exluded from sample %s because it doesnt have anisotropy correction, and other specimens are very anistropic\n"%(specimen,sample)
                                WARNING_tmp=WARNING_tmp+"excluding specimen %s; "%(specimen)
                                del tmp_Grade_A_samples[sample][specimen]
                                
                         # check if new sample pass criteria:
                        #print warning_messeage
                        if len(tmp_Grade_A_samples[sample].keys())>=self.accept_new_parameters['sample_int_n']:
                            Best_interpretations,best_mean,best_std=find_sample_min_std(tmp_Grade_A_samples[sample])
                            sample_acceptable_min,sample_acceptable_max = find_sample_min_max_interpretation (tmp_Grade_A_samples[sample],self.accept_new_parameters)
                            sample_int_interval_uT=sample_acceptable_max-sample_acceptable_min
                            sample_int_interval_perc=100*((sample_acceptable_max-sample_acceptable_min)/best_mean)       

                            # check if interpretation pass criteria (if yes ignore the specimens). if no, keep it the old way:
                            if ( self.accept_new_parameters['sample_int_sigma_uT'] ==0 and self.accept_new_parameters['sample_int_sigma_perc']==0 ) or \
                               (best_std <= self.accept_new_parameters['sample_int_sigma_uT'] or 100*(best_std/best_mean) <= self.accept_new_parameters['sample_int_sigma_perc']):
                                if sample_int_interval_uT <= self.accept_new_parameters['sample_int_interval_uT'] or sample_int_interval_perc <= self.accept_new_parameters['sample_int_interval_perc']:
                                    Grade_A_samples[sample]=copy.deepcopy(tmp_Grade_A_samples[sample])
                                    WARNING=WARNING_tmp
                                    thellier_interpreter_log.write(warning_messeage)

                                
            #--------------------------------------------------------------
            # check for outlier specimens
            # Outlier check is done only if
            # (1) number of specimen >= accept_new_parameters['sample_int_n_outlier_check']
            # (2) an outlier exists if one (and only one!) specimen has an outlier result defined
            # by:
            # Bmax(specimen_1) < mean[max(specimen_2),max(specimen_3),max(specimen_3)...] - 2*sigma
            # or
            # Bmin(specimen_1) < mean[min(specimen_2),min(specimen_3),min(specimen_3)...] + 2*sigma
            # (3) 2*sigma > 5 microT
            #--------------------------------------------------------------


            WARNING=""
            # check for outlier specimen
            exclude_specimen=""
            exclude_specimens_list=[]
            if len(Grade_A_samples[sample].keys())>=float(self.accept_new_parameters['sample_int_n_outlier_check']):
                thellier_interpreter_log.write( "-I- check outlier for sample %s \n"%sample)
                all_specimens=Grade_A_samples[sample].keys()
                for specimen in all_specimens:
                    B_min_array,B_max_array=[],[]
                    for specimen_b in all_specimens:
                        if specimen_b==specimen: continue
                        B_min_array.append(min(Grade_A_samples[sample][specimen_b]))
                        B_max_array.append(max(Grade_A_samples[sample][specimen_b]))
                    if max(Grade_A_samples[sample][specimen]) < (mean(B_min_array) - 2*std(B_min_array,ddof=1)):# and 2*std(B_min_array,ddof=1) >3.:
                        if specimen not in exclude_specimens_list:
                            exclude_specimens_list.append(specimen)
                    if min(Grade_A_samples[sample][specimen]) > (mean(B_max_array) + 2*std(B_max_array,ddof=1)):# and 2*std(B_max_array,ddof=1) >3 :
                           if specimen not in exclude_specimens_list:
                            exclude_specimens_list.append(specimen)
                         
                if len(exclude_specimens_list)>1:
                    thellier_interpreter_log.write( "-I- specimen %s outlier check: more than one specimen can be outlier. first ones are : %s,%s... \n" %(sample,exclude_specimens_list[0],exclude_specimens_list[1]))
                    exclude_specimens_list=[]

                if len(exclude_specimens_list)==1 :
                    #print exclude_specimens_list
                    exclude_specimen=exclude_specimens_list[0]
                    del Grade_A_samples[sample][exclude_specimen]
                    thellier_interpreter_log.write( "-W- WARNING: specimen %s is exluded from sample %s because of an outlier result.\n"%(exclude_specimens_list[0],sample))
                    WARNING=WARNING+"excluding specimen %s; "%(exclude_specimens_list[0])





                                        
            
            #--------------------------------------------------------------
            #  diaplay the interpretation (not only the interpretations that pass samples criteria) after the nterpreter end running.
            #--------------------------------------------------------------

            # if only one specimen take the interpretation with maximum frac
            if len(Grade_A_samples[sample].keys()) == 1:
                specimen=Grade_A_samples[sample].keys()[0]
                frac_max=0
                for TEMP in All_grade_A_Recs[specimen].keys():
                    if All_grade_A_Recs[specimen][TEMP]['specimen_frac']>frac_max:
                        best_intensity=All_grade_A_Recs[specimen][TEMP]['specimen_int_uT']
                for TEMP in All_grade_A_Recs[specimen].keys():                        
                    if All_grade_A_Recs[specimen][TEMP]['specimen_int_uT']==best_intensity:
                        self.Data[specimen]['pars'].update(All_grade_A_Recs[specimen][TEMP])
                        self.Data[specimen]['pars']['saved']=True
                        if sample not in self.Data_samples.keys():
                          self.Data_samples[sample]={}
                        self.Data_samples[sample][specimen]=self.Data[specimen]['pars']['specimen_int_uT']

            if len(Grade_A_samples[sample].keys()) > 1:
                Best_interpretations,best_mean,best_std=find_sample_min_std(Grade_A_samples[sample])
                for specimen in Grade_A_samples[sample].keys():
                    for TEMP in All_grade_A_Recs[specimen].keys():
                        if All_grade_A_Recs[specimen][TEMP]['specimen_int_uT']==Best_interpretations[specimen]:
                            self.Data[specimen]['pars'].update(All_grade_A_Recs[specimen][TEMP])
                            self.Data[specimen]['pars']['saved']=True
                            if sample not in self.Data_samples.keys():
                              self.Data_samples[sample]={}
                            self.Data_samples[sample][specimen]=self.Data[specimen]['pars']['specimen_int_uT']

        #--------------------------------------------------------------
        # calcuate STDEV-OPT 'best means' and write results to files
        #--------------------------------------------------------------

            if self.accept_new_parameters['sample_int_stdev_opt']:
                n_no_atrm=0
                
                for specimen in Grade_A_samples[sample].keys():
                    if "AniSpec" not in self.Data[specimen].keys():
                        n_no_atrm+=1
                        
                no_cooling_rate=True
                for specimen in Grade_A_samples[sample].keys():
                    if "cooling_rate_data" in self.Data[specimen].keys():
                        if "CR_correction_factor" in self.Data[specimen]["cooling_rate_data"].keys():
                            if self.Data[specimen]["cooling_rate_data"]["CR_correction_factor"]!= -1 and self.Data[specimen]["cooling_rate_data"]["CR_correction_factor"]!= -999:
                                no_cooling_rate=False
                                
     
                if len(Grade_A_samples[sample].keys())>=self.accept_new_parameters['sample_int_n']:
                    Best_interpretations,best_mean,best_std=find_sample_min_std(Grade_A_samples[sample])
                    sample_acceptable_min,sample_acceptable_max = find_sample_min_max_interpretation (Grade_A_samples[sample],self.accept_new_parameters)
                    sample_int_interval_uT=sample_acceptable_max-sample_acceptable_min
                    sample_int_interval_perc=100*((sample_acceptable_max-sample_acceptable_min)/best_mean)       
                    TEXT= "-I- sample %s 'STDEV-OPT interpretation: "%sample
                    for ss in Best_interpretations.keys():
                        TEXT=TEXT+"%s=%.1f, "%(ss,Best_interpretations[ss])
                    thellier_interpreter_log.write(TEXT+"\n")
                    thellier_interpreter_log.write("-I- sample %s STDEV-OPT mean=%f, STDEV-OPT std=%f \n"%(sample,best_mean,best_std))
                    thellier_interpreter_log.write("-I- sample %s STDEV-OPT minimum/maximum accepted interpretation  %.2f,%.2f\n" %(sample,sample_acceptable_min,sample_acceptable_max))


                    # check if interpretation pass criteria:
                    if ( self.accept_new_parameters['sample_int_sigma_uT'] ==0 and self.accept_new_parameters['sample_int_sigma_perc']==0 ) or \
                       (best_std <= self.accept_new_parameters['sample_int_sigma_uT'] or 100*(best_std/best_mean) <= self.accept_new_parameters['sample_int_sigma_perc']):
                        if sample_int_interval_uT <= self.accept_new_parameters['sample_int_interval_uT'] or sample_int_interval_perc <= self.accept_new_parameters['sample_int_interval_perc']:
                            # write the interpretation to a redo file
                            for specimen in Grade_A_samples[sample].keys():
                                #print Redo_data_specimens[specimen]
                                for TEMP in All_grade_A_Recs[specimen].keys():
                                    if All_grade_A_Recs[specimen][TEMP]['specimen_int_uT']==Best_interpretations[specimen]:
                                        t_min=All_grade_A_Recs[specimen][TEMP]['measurement_step_min']
                                        t_max=All_grade_A_Recs[specimen][TEMP]['measurement_step_max']
                                        
                                            
                                        Fout_STDEV_OPT_redo.write("%s\t%i\t%i\n"%(specimen,t_min,t_max))

                                    # write the interpretation to the specimen file
                                        #B_lab=float(Redo_data_specimens[specimen][Best_interpretations[specimen]]['lab_dc_field'])*1e6
                                        B_lab=float(All_grade_A_Recs[specimen][TEMP]['lab_dc_field'])*1e6
                                        sample=All_grade_A_Recs[specimen][TEMP]['er_sample_name']
                                        if 'AC_specimen_correction_factor' in All_grade_A_Recs[specimen][TEMP].keys():
                                            Anisotropy_correction_factor="%.2f"%float(All_grade_A_Recs[specimen][TEMP]['AC_specimen_correction_factor'])
                                        else:
                                            Anisotropy_correction_factor="-"                
                                        if  All_grade_A_Recs[specimen][TEMP]["NLT_specimen_correction_factor"] != -1:
                                            NLT_correction_factor="%.2f"%float(All_grade_A_Recs[specimen][TEMP]['NLT_specimen_correction_factor'])
                                        else:
                                            NLT_correction_factor="-"

                                        if  All_grade_A_Recs[specimen][TEMP]["specimen_int_corr_cooling_rate"] != -999 and All_grade_A_Recs[specimen][TEMP]["specimen_int_corr_cooling_rate"] != -1:
                                            CR_correction_factor="%.2f"%float(All_grade_A_Recs[specimen][TEMP]['specimen_int_corr_cooling_rate'])
                                        else:
                                            CR_correction_factor="-"

                                        Fout_STDEV_OPT_specimens.write("%s\t%s\t%.2f\t%i\t%i\t%.0f\t%s\t%s\t%s\t"\
                                                             %(sample,specimen,float(Best_interpretations[specimen]),t_min-273,t_max-273,B_lab,Anisotropy_correction_factor,NLT_correction_factor,CR_correction_factor))
                                        String=""
                                        for key in accept_specimen_keys:
                                            String=String+"%.2f"%(All_grade_A_Recs[specimen][TEMP][key])+"\t"
                                        Fout_STDEV_OPT_specimens.write(String[:-1]+"\n")
                                                     
                            # write the interpretation to the sample file
                          
                            if n_no_atrm>0:
                                 WARNING=WARNING+"% i specimens with no anisotropy correction; "%int(n_no_atrm)
                            if no_cooling_rate:
                                 WARNING=WARNING+" No cooling rate corrections; "
                            
                            String="%s\t%i\t%.2f\t%.2f\t%.3f\t%.2f\t%.2f\t%s\n"%(sample,len(Best_interpretations),best_mean,best_std,100*(best_std/best_mean),sample_int_interval_uT,sample_int_interval_perc,WARNING)
                            Fout_STDEV_OPT_samples.write(String)
                        else:
                         thellier_interpreter_log.write("-I- sample %s FAIL on sample_int_interval_uT or sample_int_interval_perc\n"%sample)                    
                    else:
                         thellier_interpreter_log.write("-I- sample %s FAIL on sample_int_sigma_uT or sample_int_sigma_perc\n"%sample)
                                                                            
        #--------------------------------------------------------------
        # calcuate Bootstarp and write results to files
        #--------------------------------------------------------------

            if self.accept_new_parameters['sample_int_bs'] or self.accept_new_parameters['sample_int_bs_par']:
               BOOTSTRAP_N=self.preferences['BOOTSTRAP_N']
               Grade_A_samples_BS={} 
               if len(Grade_A_samples[sample].keys()) >= self.accept_new_parameters['sample_int_n']:
                   for specimen in Grade_A_samples[sample].keys():
                        if specimen not in Grade_A_samples_BS.keys() and len(Grade_A_samples[sample][specimen])>0:
                           Grade_A_samples_BS[specimen]=[]
                        for B in Grade_A_samples[sample][specimen]:
                           Grade_A_samples_BS[specimen].append(B)
                        Grade_A_samples_BS[specimen].sort()
                        specimen_int_max_slope_diff=max(Grade_A_samples_BS[specimen])/min(Grade_A_samples_BS[specimen])
                        if specimen_int_max_slope_diff>self.accept_new_parameters['specimen_int_max_slope_diff']:
                           thellier_interpreter_log.write( "-I- specimen %s Failed specimen_int_max_slope_diff\n"%specimen,Grade_A_samples_BS[specimen])
                           del Grade_A_samples_BS[specimen]
                
               if len(Grade_A_samples_BS.keys())>=self.accept_new_parameters['sample_int_n']:
        
                   BS_means_collection=[]
                   for i in range(BOOTSTRAP_N):
                       B_BS=[]
                       for j in range(len(Grade_A_samples_BS.keys())):
                           LIST=list(Grade_A_samples_BS.keys())
                           specimen=random.choice(LIST)
                           if self.accept_new_parameters['sample_int_bs']:
                               B=random.choice(Grade_A_samples_BS[specimen])
                           if self.accept_new_parameters['sample_int_bs_par']:
                               B=random.uniform(min(Grade_A_samples_BS[specimen]),max(Grade_A_samples_BS[specimen]))
                           B_BS.append(B)
                       BS_means_collection.append(mean(B_BS))
                       
                   BS_means=array(BS_means_collection)
                   BS_means.sort()
                   sample_median=median(BS_means)
                   sample_std=std(BS_means,ddof=1)
                   sample_68=[BS_means[(0.16)*len(BS_means)],BS_means[(0.84)*len(BS_means)]]
                   sample_95=[BS_means[(0.025)*len(BS_means)],BS_means[(0.975)*len(BS_means)]]


                   thellier_interpreter_log.write( "-I-  bootstrap mean sample %s: median=%f, std=%f\n"%(sample,sample_median,sample_std))
                   String="%s\t%i\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%s\n"%\
                           (sample,len(Grade_A_samples[sample].keys()),sample_median,sample_68[0],sample_68[1],sample_95[0],sample_95[1],sample_std,100*(sample_std/sample_median),WARNING)
                   if self.accept_new_parameters['sample_int_bs']:
                       Fout_BS_samples.write(String)
                   if self.accept_new_parameters['sample_int_bs_par']:
                       Fout_BS_PAR_samples.write(String)


                                                  
        
        thellier_interpreter_log.write( "-I- Statistics:\n")
        thellier_interpreter_log.write( "-I- number of specimens analzyed = %i\n" % len(specimens_list)  )
        thellier_interpreter_log.write( "-I- number of sucsessful 'acceptable' specimens = %i\n" % len(All_grade_A_Recs.keys()))   

        runtime_sec = time.time() - start_time
        m, s = divmod(runtime_sec, 60)
        h, m = divmod(m, 60)
        thellier_interpreter_log.write( "-I- runtime hh:mm:ss is " + "%d:%02d:%02d\n" % (h, m, s))
        thellier_interpreter_log.write( "-I- Finished sucsessfuly.\n")
        thellier_interpreter_log.write( "-I- DONE\n")


        # close all files

        thellier_interpreter_log.close()
        thellier_interpreter_all.close()
        Fout_specimens_bounds.close()
        if self.accept_new_parameters['sample_int_stdev_opt']: 
            Fout_STDEV_OPT_redo.close()
            Fout_STDEV_OPT_specimens.close()
            Fout_STDEV_OPT_samples.close()
        if self.accept_new_parameters['sample_int_bs']:
            Fout_BS_samples.close()
        if self.accept_new_parameters['sample_int_bs_par']:
            Fout_BS_PAR_samples.close()
            
        os.system('\a')
        dlg1 = wx.MessageDialog(self,caption="Message:", message="Interpreter finished sucsessfuly\nCheck output files in folder /thellier_interpreter in the current project directory" ,style=wx.OK|wx.ICON_INFORMATION)


        # display the interpretation of the current specimen:
        self.pars=self.Data[self.s]['pars']
        self.clear_boxes()
        self.draw_figure(self.s)
        self.update_GUI_with_new_interpretation()


        dlg1.ShowModal()
        dlg1.Destroy()
        busy_frame.Destroy()
        """
    #----------------------------------------------------------------------
    # END OF  calculate anisotropy_tensors()

    #----------------------------------------------------------------------
        
    def read_redo_file(self,redo_file):
        """
        Read previous interpretation from a redo file
        and update gui with the new interpretation
        """
        print ("-I- read redo file and processing new temperature bounds")
        not_called = """print "calling read_redo_file()"
        self.redo_specimens={}
        # first delete all previous interpretation
        for sp in self.Data.keys():
            del self.Data[sp]['pars']
            self.Data[sp]['pars']={}
            self.Data[sp]['pars']['lab_dc_field']=self.Data[sp]['lab_dc_field']
            self.Data[sp]['pars']['er_specimen_name']=self.Data[sp]['er_specimen_name']   
            self.Data[sp]['pars']['er_sample_name']=self.Data[sp]['er_sample_name']
            #print sp
            #print self.Data[sp]['pars']
        self.Data_samples={}
        
        fin=open(redo_file,'rU')
        for Line in fin.readlines():
          line=Line.strip('\n').split()
          specimen=line[0]
          tmin_kelvin=float(line[1])
          tmax_kelvin=float(line[2])
          if specimen not in self.redo_specimens.keys():
            self.redo_specimens[specimen]={}
          self.redo_specimens[specimen]['t_min']=float(tmin_kelvin)
          self.redo_specimens[specimen]['t_max']=float(tmax_kelvin)
          if specimen in self.Data.keys():
              if tmin_kelvin not in self.Data[specimen]['t_Arai'] or tmax_kelvin not in self.Data[specimen]['t_Arai'] :
                  print ("-W- WARNING: cant fit temperature bounds in the redo file to the actual measurement. specimen %s\n"%specimen)
              else:
                  try:
                      self.Data[specimen]['pars']=self.get_PI_parameters(specimen,float(tmin_kelvin),float(tmax_kelvin))
                      self.Data[specimen]['pars']['saved']=True
                      # write intrepretation into sample data
                      sample=self.Data_hierarchy['specimens'][specimen]
                      if sample not in self.Data_samples.keys():
                          self.Data_samples[sample]={}
                      self.Data_samples[sample][specimen]=self.Data[specimen]['pars']['specimen_int_uT']
                  except:
                      print ("-E- ERROR. Cant calculate PI paremeters for specimen %s using redo file. Check!"%(specimen))
          else:
              print ("-W- WARNING: Cant find specimen %s from redo file in measurement file!\n"%specimen)
        fin.close()
        self.pars=self.Data[self.s]['pars']
        self.clear_boxes()
        self.draw_figure(self.s)
        self.update_GUI_with_new_interpretation()"""

    #----------------------------------------------------------------------            

    def write_acceptance_criteria_to_file(self):
        print "calling write_acceptance_criteria_to_file()"
 #       import copy
        """
        Write new acceptance criteria to pmag_criteria.txt
        """
        # check if an old pmag_criteria.txt exist:
        not_called = """other_criteria={}
        try:
            fin=open(self.WD+"/"+"pmag_criteria.txt",'rU')
            lines=""
            line=fin.readline()
            line=fin.readline()
            header=line.strip('\n').split()
            code_index=header.index("pmag_criteria_code")

            for line in fin.readlines():
                code=line[code_index]
                if "IE-" not in code:
                    for i in range(len(header)):
                        if line[i]!="":
                            try:
                                float(line[i])
                            except:
                                continue
                            other_criteria[code][header[i]]=float(line[i])
        except:
             pass
    
                
            
        fout=open(self.WD+"/"+"pmag_criteria.txt",'w')
        String="tab\tpmag_criteria\n"
        fout.write(String)
        specimen_criteria_list=self.criteria_list+["specimen_int_max_slope_diff"]+['check_aniso_ftest']+['anisotropy_alt']
        sample_criteria_list=[key for key in self.accept_new_parameters.keys() if "sample" in key]
        if self.accept_new_parameters['sample_int_stdev_opt'] == True:                                      
            for k in ['sample_int_bs','sample_int_bs_par','sample_int_BS_68_uT','sample_int_BS_68_perc','sample_int_BS_95_uT','sample_int_BS_95_perc']:
                sample_criteria_list.remove(k)
                if "specimen_int_max_slope_diff" in specimen_criteria_list:
                    specimen_criteria_list.remove("specimen_int_max_slope_diff")

        else:
            for k in ['sample_int_sigma_uT','sample_int_sigma_perc','sample_int_interval_uT','sample_int_stdev_opt','sample_aniso_threshold_perc']:
                sample_criteria_list.remove(k)
        for k in ['sample_int_sigma_uT','sample_int_sigma_perc','sample_int_interval_uT','sample_int_interval_perc','sample_aniso_threshold_perc','sample_int_BS_68_uT','sample_int_BS_68_perc','sample_int_BS_95_uT','sample_int_BS_95_perc',]:
            if  k in sample_criteria_list:
                if float(self.accept_new_parameters[k]) > 999:
                    sample_criteria_list.remove(k)
        if  float(self.accept_new_parameters["sample_int_n_outlier_check"])> 99:
            sample_criteria_list.remove("sample_int_n_outlier_check")

        if "sample_int_sigma_uT"  in sample_criteria_list and "sample_int_sigma" not in sample_criteria_list:
            sample_criteria_list.append("sample_int_sigma")
            self.accept_new_parameters["sample_int_sigma"]=float(self.accept_new_parameters["sample_int_sigma_uT"])*1e-6
        
        if "specimen_int_max_slope_diff" in  specimen_criteria_list:
            if float(self.accept_new_parameters['specimen_int_max_slope_diff'])>999:
                specimen_criteria_list.remove("specimen_int_max_slope_diff")
        c_list=copy.copy(specimen_criteria_list)       
        for criteria in c_list:
            if criteria in (self.high_threshold_velue_list + ['anisotropy_alt']) and float(self.accept_new_parameters[criteria])>100:
                specimen_criteria_list.remove(criteria)
            #if criteria in ['specimen_g'] and float(self.accept_new_parameters[criteria])>100:
            if criteria in self.low_threshold_velue_list and float(self.accept_new_parameters[criteria])<0.1:
                specimen_criteria_list.remove(criteria)                

        header="pmag_criteria_code\t"
        for key in sample_criteria_list:
            header=header+key+"\t"
        for key in specimen_criteria_list:
            header=header+key+"\t"
        header=header+"specimen_scat\t"

        # other criteria (not paleointensity)
        for code in other_criteria.keys():
            for key in other_criteria[code].keys():
                header=header+key+"\t"
        fout.write(header[:-1]+"\n")
                    
        line="IE-SPEC:IE-SAMP\t"
        for key in sample_criteria_list:
            if key in ['sample_int_bs','sample_int_bs_par','sample_int_stdev_opt','check_aniso_ftest']:
                line=line+"%s"%str(self.accept_new_parameters[key])+"\t"
            elif key in ['sample_int_sigma']:
                line=line+"%.2e"%self.accept_new_parameters[key]+"\t"                
            else:
                line=line+"%f"%self.accept_new_parameters[key]+"\t"
                
        for key in specimen_criteria_list:
            if key=='check_aniso_ftest':
                line=line+str(self.accept_new_parameters[key])+"\t"
            else:    
                line=line+"%f"%self.accept_new_parameters[key]+"\t"
        if self.accept_new_parameters["specimen_scat"]:
            line=line+"True"+"\t"
        else:
            line=line+"False"+"\t"

        # other criteria (not paleointensity)
        for code in other_criteria.keys():
            for key in other_criteria[code].keys():
                line=line+other_criteria[code][key]+"\t"
    
        fout.write(line[:-1]+"\n")
        fout.close()
            
    #----------------------------------------------------------------------            
      """
    

    def on_menu_results_data (self, event):

        # Results of all the samples that passed the criteria
        
        # search for ages and Latitudes
        print "calling on_menu_results_data()"
        not_called="""
        samples_list=self.Data_samples.keys()
        samples_list.sort()
        Results_table_data={}
        for sample in samples_list:

            Age,age_unit,age_range_low,age_range_high="","","",""
            lat,lon,VADM,VADM_sigma="","","",""

            found_age,found_lat=False,False

            # Find the mean paleointenisty for each sample
            tmp_B=[]
            for spec in self.Data_samples[sample].keys():
                tmp_B.append( self.Data_samples[sample][spec])
            if len(tmp_B)<1:
                continue
            tmp_B=array(tmp_B)
            B_uT=mean(tmp_B)
            B_std_uT=std(tmp_B,ddof=1)
            B_std_perc=100*(B_std_uT/B_uT)
            
            # check if sample passed the criteria
            sample_pass_criteria=False
            if len(tmp_B)>=self.accept_new_parameters['sample_int_n']:
                if (self.accept_new_parameters['sample_int_sigma_uT']==0 and self.accept_new_parameters['sample_int_sigma_perc']==0) or\
                   ( B_std_uT <=self.accept_new_parameters['sample_int_sigma_uT'] or B_std_perc <= self.accept_new_parameters['sample_int_sigma_perc']):
                    if ( (max(tmp_B)-min(tmp_B)) <= self.accept_new_parameters['sample_int_interval_uT'] or 100*((max(tmp_B)-min(tmp_B))/mean((tmp_B))) <= self.accept_new_parameters['sample_int_interval_perc']):
                        sample_pass_criteria=True

            if not sample_pass_criteria:
                continue

            Results_table_data[sample]={}
            
            # search for samples age in er_ages.txt by sample or by site
            site = self.Data_info["er_samples"][sample]['er_site_name']
            found_age=False
            if sample in self.Data_info["er_ages"].keys():
                age_key=sample
            elif site in self.Data_info["er_ages"].keys():
                age_key=site
            else:
                age_key=""
            if age_key !="":
                try:
                    age_unit=self.Data_info["er_ages"][age_key]["age_unit"]                
                except:
                    age_unit="unknown"               
                    
                if self.Data_info["er_ages"][age_key]["age"] !="":
                    Age = float(self.Data_info["er_ages"][age_key]["age"])
                    found_age=True
                    
                if "age_range_low" in self.Data_info["er_ages"][age_key].keys() and "age_range_high" in self.Data_info["er_ages"][age_key].keys():
                   age_range_low=float(self.Data_info["er_ages"][age_key]["age_range_low"])
                   age_range_high=float(self.Data_info["er_ages"][age_key]["age_range_high"])
                   
                   if not found_age:
                       Age=(age_range_low+age_range_high)/2
                       found_age=True

                elif "age_sigma" in self.Data_info["er_ages"][age_key].keys() and found_age:
                   age_range_low=Age-float(self.Data_info["er_ages"][age_key]["age_sigma"])
                   age_range_high= Age+float(self.Data_info["er_ages"][age_key]["age_sigma"])

                elif found_age:
                   age_range_low=Age
                   age_range_high=Age

            # convert ages from Years BP to Years Cal AD (+/-)
                if "Years BP" in age_unit:
                    Age=1950-Age
                    age_range_low=1950-age_range_low
                    age_range_high=1950-age_range_high
                    age_unit="Years Cal AD (+/-)"             
            
            # search for Lon/Lat
            if sample in self.Data_info["er_samples"].keys() and "site_lat" in self.Data_info["er_samples"][sample].keys():
                lat=float(self.Data_info["er_samples"][sample]["site_lat"])
                lon=float(self.Data_info["er_samples"][sample]["site_lon"])
                found_lat=True
                
            elif site in self.Data_info["er_sites"].keys() and "site_lat" in self.Data_info["er_sites"][site].keys():
                lat=float(self.Data_info["er_sites"][site]["site_lat"])
                lon=float(self.Data_info["er_sites"][site]["site_lon"])
                found_lat=True

            if found_lat:
                VADM=self.b_vdm(B_uT*1e-6,lat)*1e-21
                VADM_plus=self.b_vdm((B_uT+B_std_uT)*1e-6,lat)*1e-21
                VADM_minus=self.b_vdm((B_uT-B_std_uT)*1e-6,lat)*1e-21
                VADM_sigma=(VADM_plus-VADM_minus)/2
                
            Results_table_data[sample]["N"]="%i"%(len(tmp_B))            
            Results_table_data[sample]["B_uT"]="%.1f"%(B_uT)
            Results_table_data[sample]["B_std_uT"]="%.1f"%(B_std_uT)
            Results_table_data[sample]["B_std_perc"]="%.1f"%(B_std_perc)
            if found_lat:
                Results_table_data[sample]["Lat"]="%f"%lat
                Results_table_data[sample]["Lon"]="%f"%lon
                Results_table_data[sample]["VADM"]="%.1f"%VADM
                Results_table_data[sample]["VADM_sigma"]="%.1f"%VADM_sigma
            else:
                Results_table_data[sample]["Lat"]=""
                Results_table_data[sample]["Lon"]=""
                Results_table_data[sample]["VADM"]=""
                Results_table_data[sample]["VADM_sigma"]=""
            if found_age:
                Results_table_data[sample]["Age"]="%.2f"%Age
                Results_table_data[sample]["Age_low"]="%.2f"%age_range_low
                Results_table_data[sample]["Age_high"]="%.2f"%age_range_high
            else:
                Results_table_data[sample]["Age"]=""
                Results_table_data[sample]["Age_low"]=""
                Results_table_data[sample]["Age_high"]=""
            Results_table_data[sample]["Age_units"]=age_unit
                
        samples_list= Results_table_data.keys()
        samples_list.sort()
        if len(samples_list) <1:
            return
        else:
            fout=open(self.WD+"/results_table.txt",'w')
            Keys=["sample","Lat","Lon","Age","Age_low","Age_high","Age_units","N","B_uT","B_std_uT","VADM","VADM_sigma"]
            fout.write("\t".join(Keys)+"\n")
            for sample in samples_list:
                String=sample+"\t"
                for k in Keys[1:]:
                    String=String+Results_table_data[sample][k]+"\t"
                fout.write(String[:-1]+"\n")
            fout.close()

            dlg1 = wx.MessageDialog(self,caption="Message:", message="Output results table is saved in 'results_table.txt'" ,style=wx.OK|wx.ICON_INFORMATION)
            dlg1.ShowModal()
            dlg1.Destroy()
            
        return
    """
    #----------------------------------------------------------------------            
    
        

    def read_magic_model (self):
        # Read MagIC Data model:
        print "calling read_magic_model()"
        not_called="""
        self.MagIC_model={}
        self.MagIC_model["specimens"]={}
        self.MagIC_model["er_samples"]={}
        self.MagIC_model["er_sites"]={}
        self.MagIC_model["er_locations"]={}
        self.MagIC_model["er_ages"]={}
        fail=[]
        self.MagIC_model["specimens"]=self.read_magic_file(self.WD+"/er_specimens.txt",1,'er_specimen_name')
        try:
            self.MagIC_model["specimens"]=self.read_magic_file(self.WD+"/er_specimens.txt",1,'er_specimen_name')
        except:
            print ("-W- Cant find er_specimens.txt in project directory")
            fail.append("er_specimens.txt")
            pass
        try:
            self.MagIC_model["er_samples"]=self.read_magic_file(self.WD+"/er_samples.txt",1,'er_sample_name')
        except:
            print ("-W- Cant find er_sample.txt in project directory")
            fail.append("er_sample.txt")
            pass
        try:
            self.MagIC_model["er_sites"]=self.read_magic_file(self.WD+"/er_sites.txt",1,'er_site_name')
        except:
            print ("-W- Cant find er_sites.txt in project directory")
            fail.append("er_sites.txt")
            pass
        try:
            self.MagIC_model["er_locations"]=self.read_magic_file(self.WD+"/er_locations.txt",1,'er_location_name')
        except:
            print ("-W- Cant find er_locations.txt in project directory")
            fail.append("er_locations.txt")
            pass

        try:
            self.MagIC_model["er_ages"]=self.read_magic_file(self.WD+"/er_ages",1,'er_site_name')
        except:
            print ("-W- Cant find er_ages.txt in project directory")
            pass

        return (fail)"""

                          
    def classy_read_magic_file(self,path,ignore_lines_n,sort_by_this_name): # only called for 'pmag_specimens.txt'
        print "calling read_magic_file() in thellier_gui_spd_lj.py"
        print path
        DATA={}
        fin=open(path,'rU')
        #ignore first lines
        for i in range(ignore_lines_n):
            fin.readline()
        #header
        line=fin.readline()
        header=line.strip('\n').split('\t')
        #print header
        for line in fin.readlines():
            if line[0]=="#":
                continue
            else: print "line[0] != '#'"
            tmp_data={}
            tmp_line=line.strip('\n').split('\t')
            #print tmp_line
            for i in range(len(tmp_line)):
                if i>= len(header):
                    continue
                else: print "something tripped"
                tmp_data[header[i]]=tmp_line[i]
            DATA[tmp_data[sort_by_this_name]]=tmp_data
        fin.close()        
        print "Data from read_magic_file in nothing:  ", DATA
        return(DATA)


#===========================================================
# calculate PI statistics
#===========================================================



    def get_new_T_PI_parameters(self,event):
        
        """
        calcualte statisics when temperatures are selected
        """
        print "calling get_new_T_PI_parameters"
        #remember the last saved interpretation
        not_called = """
        if "saved" in self.pars.keys():
            if self.pars['saved']:
                self.last_saved_pars={}
                for key in self.pars.keys():
                    self.last_saved_pars[key]=self.pars[key]
        self.pars['saved']=False
        t1=self.tmin_box.GetStringSelection()
        t2=self.tmax_box.GetStringSelection()

        if (t1 == "" or t2==""):
          return()
        if float(t2) < float(t1):
          return()

        if float(t2) < float(t1):
          return()

        index_1=self.T_list.index(t1)
        index_2=self.T_list.index(t2)

       
        if (index_2-index_1)+1 >= self.accept_new_parameters['specimen_int_n']:
            if self.Data[self.s]['T_or_MW']!="MW":
                self.pars=self.get_PI_parameters(self.s,float(t1)+273.,float(t2)+273.)
            else:
                self.pars=self.get_PI_parameters(self.s,float(t1),float(t2))
                
            self.update_GUI_with_new_interpretation()
      """
    

    def get_PI_parameters(self,s,tmin,tmax):
        print "calling get_PI_parameters() from thellier_gui_spd_lj.py"
        print "self", self, str(self.Data)[:500] + "..."


        def calculate_ftest(s,sigma,nf):
            print "calling calculate_ftest() in get_PI_parameters()"
            chibar=(s[0][0]+s[1][1]+s[2][2])/3.
            t=array(linalg.eigvals(s))
            F=0.4*(t[0]**2+t[1]**2+t[2]**2 - 3*chibar**2)/(float(sigma)**2)

            return(F)

        """
        calcualte statisics 
        """
        
        print "calculating statistics"
        pars=self.Data[s]['pars']
        datablock = self.Data[s]['datablock']
        pars=self.Data[s]['pars']
        print "new pars", pars
        # get MagIC mothod codes:

        #pars['magic_method_codes']="LP-PI-TRM" # thellier Method
        
        
        t_Arai=self.Data[s]['t_Arai']
        x_Arai=self.Data[s]['x_Arai']
        y_Arai=self.Data[s]['y_Arai']
        x_tail_check=self.Data[s]['x_tail_check']
        y_tail_check=self.Data[s]['y_tail_check']

        zijdblock=self.Data[s]['zijdblock']        
        z_temperatures=self.Data[s]['z_temp']
        print "got through a few stats"
        #print tmin,tmax,z_temperatures
        # check tmin
        if tmin not in t_Arai or tmin not in z_temperatures:
            print "t_Arai   ", t_Arai
            print "z_temperatures  ", z_temperatures
            print "tmin not in t_Arai or not in z_temperatures"
            return(pars)
        
        # check tmax
        if tmax not in t_Arai or tmin not in z_temperatures:
            print "returning pars because tmax not in t_Arai"
            return(pars)

        print "got past first 2 if statements"

        start=t_Arai.index(tmin)
        end=t_Arai.index(tmax)

        if end-start < float(self.accept_new_parameters['specimen_int_n'] -1):
          print "returning pars because??"
          return(pars)
                                                 
        #-------------------------------------------------
        # calualte PCA of the zerofield steps
        # MAD calculation following Kirschvink (1980)
        # DANG following Tauxe and Staudigel (2004)
        #-------------------------------------------------               
         
        pars["measurement_step_min"]=float(tmin)
        pars["measurement_step_max"]=float(tmax)
 
        zstart=z_temperatures.index(tmin)
        zend=z_temperatures.index(tmax)

        zdata_segment=self.Data[s]['zdata'][zstart:zend+1]

        print "elephant"
        #  PCA in 2 lines
        M = (zdata_segment-mean(zdata_segment.T,axis=1)).T # subtract the mean (along columns)
        [eigenvalues,eigenvectors] = linalg.eig(cov(M)) # attention:not always sorted

        # sort eigenvectors and eigenvalues
        eigenvalues=list(eigenvalues)
        tmp=[0,1,2]
        t1=max(eigenvalues);index_t1=eigenvalues.index(t1);tmp.remove(index_t1)
        t3=min(eigenvalues);index_t3=eigenvalues.index(t3);tmp.remove(index_t3)
        index_t2=tmp[0];t2=eigenvalues[index_t2]
        v1=real(array(eigenvectors[:,index_t1]))
        v2=real(array(eigenvectors[:,index_t2]))
        v3=real(array(eigenvectors[:,index_t3]))

        # chech if v1 is the "right" polarity
        cm=array(mean(zdata_segment.T,axis=1)) # center of mass
        v1_plus=v1*sqrt(sum(cm**2))
        v1_minus=v1*-1*sqrt(sum(cm**2))
        test_v=zdata_segment[0]-zdata_segment[-1]

        if sqrt(sum((v1_minus-test_v)**2)) < sqrt(sum((v1_plus-test_v)**2)):
         DIR_PCA=self.cart2dir(v1*-1)
         best_fit_vector=v1*-1
        else:
         DIR_PCA=self.cart2dir(v1)
         best_fit_vector=v1

        # MAD Kirschvink (1980)
        MAD=math.degrees(arctan(sqrt((t2+t3)/t1)))

        # DANG Tauxe and Staudigel 2004
        DANG=math.degrees( arccos( ( dot(cm, best_fit_vector) )/( sqrt(sum(cm**2)) * sqrt(sum(best_fit_vector**2)))))


        # best fit PCA direction
        pars["specimen_dec"] =  DIR_PCA[0]
        pars["specimen_inc"] =  DIR_PCA[1]
        pars["specimen_PCA_v1"] =best_fit_vector
        if t1 <0 or t1==0:
            t1=1e-10
        if t2 <0 or t2==0:
            t2=1e-10
        if t3 <0 or t3==0:
            t3=1e-10
            
        pars["specimen_PCA_sigma_max"] =  sqrt(t1)
        pars["specimen_PCA_sigma_int"] =  sqrt(t2)
        pars["specimen_PCA_sigma_min"] =  sqrt(t3)
            

        # MAD Kirschvink (1980)
        pars["specimen_int_mad"]=MAD
        pars["specimen_dang"]=DANG


        #-------------------------------------------------
        # calualte PCA of the pTRMs over the entire temperature range
        # and calculate the angular difference to the lab field
        # MAD calculation following Kirschvink (1980)
        #-------------------------------------------------
        
        PTRMS = self.Data[s]['PTRMS'][1:]
        CART_pTRMS_orig=array([self.dir2cart(row[1:4]) for row in PTRMS])
        #CART_pTRMS=[row/sqrt(sum((array(row)**2))) for row in CART_pTRMS_orig]
##        print "CART_pTRMS_orig",CART_pTRMS_orig
##        print "----"
        
        #  PCA in 2 lines
        M = (CART_pTRMS_orig-mean(CART_pTRMS_orig.T,axis=1)).T # subtract the mean (along columns)
        [eigenvalues,eigenvectors] = linalg.eig(cov(M)) # attention:not always sorted

        # sort eigenvectors and eigenvalues
        eigenvalues=list(eigenvalues)
        tmp=[0,1,2]
        t1=max(eigenvalues);index_t1=eigenvalues.index(t1);tmp.remove(index_t1)
        t3=min(eigenvalues);index_t3=eigenvalues.index(t3);tmp.remove(index_t3)
        index_t2=tmp[0];t2=eigenvalues[index_t2]
        v1=real(array(eigenvectors[:,index_t1]))
        v2=real(array(eigenvectors[:,index_t2]))
        v3=real(array(eigenvectors[:,index_t3]))

        # chech if v1 is the "right" polarity
        cm=array(mean(CART_pTRMS_orig.T,axis=1)) # center of mass
        v1_plus=v1*sqrt(sum(cm**2))
        v1_minus=v1*-1*sqrt(sum(cm**2))
        test_v=CART_pTRMS_orig[0]-CART_pTRMS_orig[-1]

        if sqrt(sum((v1_minus-test_v)**2)) > sqrt(sum((v1_plus-test_v)**2)):
         DIR_PCA=self.cart2dir(v1*-1)
         best_fit_vector=v1*-1
        else:
         DIR_PCA=self.cart2dir(v1)
         best_fit_vector=v1

        # MAD Kirschvink (1980)
        MAD=math.degrees(arctan(sqrt((t2+t3)/t1)))


        # best fit PCA direction
        pars["specimen_ptrms_dec"] =  DIR_PCA[0]
        pars["specimen_ptrms_inc"] =  DIR_PCA[1]
        pars["specimen_ptrms_mad"]=MAD
        B_lab_unit=self.dir2cart([ self.Data[s]['Thellier_dc_field_phi'], self.Data[s]['Thellier_dc_field_theta'],1])
        pars["specimen_ptrms_angle"]=math.degrees(math.acos(dot(best_fit_vector,B_lab_unit)/(sqrt(sum(best_fit_vector**2)) * sqrt(sum(B_lab_unit**2)))))

##        print "specimen_ptrms_dec",pars["specimen_ptrms_dec"]
##        print "specimen_ptrms_inc",pars["specimen_ptrms_inc"]
##        print "B_lab_unit,v1",B_lab_unit,v1
##        print "specimen_ptrms_angle", pars["specimen_ptrms_angle"]

##        #-------------------------------------------------                     
##        # Calculate the new 'MAD box' parameter
##        # all datapoints should be inside teh M"AD box"
##        # defined by the threshold value of MAD
##        # For definitionsee Shaar and Tauxe (2012)
##        #-------------------------------------------------                     
##
##        pars["specimen_mad_scat"]="Pass"
##        self.accept_new_parameters['specimen_mad_scat']=True
##        if 'specimen_mad_scat' in self.accept_new_parameters.keys() and 'specimen_int_mad' in self.accept_new_parameters.keys() :
##            if self.accept_new_parameters['specimen_mad_scat']==True or self.accept_new_parameters['specimen_mad_scat'] in [1,"True","TRUE",'1']:
##
##                # center of mass 
##                CM_x=mean(zdata_segment[:,0])
##                CM_y=mean(zdata_segment[:,1])
##                CM_z=mean(zdata_segment[:,2])
##                CM=array([CM_x,CM_y,CM_z])
##
##                # threshold value for the distance of the point from a line:
##                # this is depends of MAD
##                # if MAD= tan-1 [ sigma_perpendicular / sigma_max ]
##                # then:
##                # sigma_perpendicular_threshold=tan(MAD_threshold)*sigma_max
##                sigma_perpendicular_threshold=abs(tan(radians(self.accept_new_parameters['specimen_int_mad'])) *  pars["specimen_PCA_sigma_max"] )
##                
##                # Line from
##                #print "++++++++++++++++++++++++++++++++++"
##                
##                for P in zdata_segment:
##                    # Find the line  P_CM that connect P to the center of mass
##                    #print "P",P
##                    #print "CM",CM
##                    P_CM=P-CM
##                    #print "P_CM",P_CM
##                    
##                    #  the dot product of vector P_CM with the unit direction vector of the best-fit liene. That's the projection of P_CM on the PCA line 
##                    best_fit_vector_unit=best_fit_vector/sqrt(sum(best_fit_vector**2))
##                    #print "best_fit_vector_unit",best_fit_vector_unit
##                    CM_P_projection_on_PCA_line=dot(best_fit_vector_unit,P_CM)
##                    #print "CM_P_projection_on_PCA_line",CM_P_projection_on_PCA_line
##
##                    # Pythagoras
##                    P_CM_length=sqrt(sum((P_CM)**2))
##                    Point_2_PCA_Distance=sqrt((P_CM_length**2-CM_P_projection_on_PCA_line**2))
##                    #print "Point_2_PCA_Distance",Point_2_PCA_Distance
##
##
##                    #print "sigma_perpendicular_threshold*2",sigma_perpendicular_threshold*2
##                    if Point_2_PCA_Distance > sigma_perpendicular_threshold*2:
##                        pars["specimen_mad_scat"]="Fail"
##                        index=999
##                        for i in range(len(self.Data[s]['zdata'])):
##                        
##                            if P[0] == self.Data[s]['zdata'][i][0] and P[1] == self.Data[s]['zdata'][i][1] and P[2] == self.Data[s]['zdata'][i][2]:
##                                index =i
##                                break
##                        #print "specimen  %s fail on mad_scat,%i"%(s,index)
##                        
##                    
##                    
##                    #CM_P_projection_on_PCA_line_length=sqrt(sum((CM_P_projection_on_PCA_line_length)**2))
        

        #-------------------------------------------------
        # York regresssion (York, 1967) following Coe (1978)
        # calculate f,fvds,
        # modified from pmag.py
        #-------------------------------------------------               

        x_Arai_segment= x_Arai[start:end+1]
        y_Arai_segment= y_Arai[start:end+1]

        x_Arai_mean=mean(x_Arai_segment)
        y_Arai_mean=mean(y_Arai_segment)

        # equations (2),(3) in Coe (1978) for b, sigma
        n=end-start+1
        x_err=x_Arai_segment-x_Arai_mean
        y_err=y_Arai_segment-y_Arai_mean

        # York b
        york_b=-1* sqrt( sum(y_err**2) / sum(x_err**2) )

        # york sigma
        york_sigma= sqrt ( (2 * sum(y_err**2) - 2*york_b*sum(x_err*y_err)) / ( (n-2) * sum(x_err**2) ) )

        # beta  parameter                
        beta_Coe=abs(york_sigma/york_b)

        # y_T is the intercept of the extrepolated line
        # through the center of mass (see figure 7 in Coe (1978))
        y_T = y_Arai_mean - york_b* x_Arai_mean

        # calculate the extarplated data points for f and fvds
        # (see figure 7 in Coe (1978))

        x_tag=(y_Arai_segment - y_T ) / york_b
        y_tag=york_b*x_Arai_segment + y_T

        # intersect of the dashed square and the horizontal dahed line  next to delta-y-5 in figure 7, Coe (1978)
        x_prime=(x_Arai_segment+x_tag) / 2
        y_prime=(y_Arai_segment+y_tag) / 2

        f_Coe=abs((y_prime[0]-y_prime[-1])/y_T)

        f_vds=abs((y_prime[0]-y_prime[-1])/self.Data[s]['vds'])

        g_Coe= 1 - (sum((y_prime[:-1]-y_prime[1:])**2) / sum((y_prime[:-1]-y_prime[1:]))**2 )

        q_Coe=abs(york_b)*f_Coe*g_Coe/york_sigma


        count_IZ= self.Data[self.s]['steps_Arai'].count('IZ')
        count_ZI= self.Data[self.s]['steps_Arai'].count('ZI')
        if count_IZ >1 and count_ZI >1:
            pars['magic_method_codes']="LP-PI-BT-IZZI"
        elif count_IZ <1 and count_ZI >1:
            pars['magic_method_codes']="LP-PI-ZI"
        elif count_IZ >1 and count_ZI <1:
            pars['magic_method_codes']="LP-PI-IZ"            
        else:
            pars['magic_method_codes']=""
            
        pars['specimen_int_n']=end-start+1
        pars["specimen_b"]=york_b
        pars["specimen_YT"]=y_T       
        pars["specimen_b_sigma"]=york_sigma
        pars["specimen_b_beta"]=beta_Coe
        pars["specimen_f"]=f_Coe
        pars["specimen_fvds"]=f_vds
        pars["specimen_g"]=g_Coe
        pars["specimen_q"]=q_Coe
        pars["specimen_int"]=-1*pars['lab_dc_field']*pars["specimen_b"]
        pars['magic_method_codes']+=":IE-TT"
        if 'x_ptrm_check' in self.Data[self.s].keys():
            if len(self.Data[self.s]['x_ptrm_check'])>0:
                pars['magic_method_codes']+=":LP-PI-ALT-PTRM"
        if 'x_tail_check' in self.Data[self.s].keys():
            if len(self.Data[self.s]['x_tail_check'])>0:
                pars['magic_method_codes']+=":LP-PI-BT-MD"


        #-------------------------------------------------
        # pTRM checks:
        # DRAT ()
        # and
        # DRATS (Tauxe and Staudigel 2004)
        #-------------------------------------------------

        x_ptrm_check_in_0_to_end,y_ptrm_check_in_0_to_end,x_Arai_compare=[],[],[]
        x_ptrm_check_in_start_to_end,y_ptrm_check_in_start_to_end=[],[]
        x_ptrm_check_for_SCAT,y_ptrm_check_for_SCAT=[],[]

        stop_scat_collect=False
        for k in range(len(self.Data[s]['ptrm_checks_temperatures'])):
          if self.Data[s]['ptrm_checks_temperatures'][k]<pars["measurement_step_max"] and self.Data[s]['ptrm_checks_temperatures'][k] in t_Arai:
            x_ptrm_check_in_0_to_end.append(self.Data[s]['x_ptrm_check'][k])
            y_ptrm_check_in_0_to_end.append(self.Data[s]['y_ptrm_check'][k])
            x_Arai_index=t_Arai.index(self.Data[s]['ptrm_checks_temperatures'][k])
            x_Arai_compare.append(x_Arai[x_Arai_index])
            if self.Data[s]['ptrm_checks_temperatures'][k]>=pars["measurement_step_min"]:
                x_ptrm_check_in_start_to_end.append(self.Data[s]['x_ptrm_check'][k])
                y_ptrm_check_in_start_to_end.append(self.Data[s]['y_ptrm_check'][k])
          if self.Data[s]['ptrm_checks_temperatures'][k] >= pars["measurement_step_min"] and self.Data[s]['ptrm_checks_starting_temperatures'][k] <= pars["measurement_step_max"] :
                x_ptrm_check_for_SCAT.append(self.Data[s]['x_ptrm_check'][k])
                y_ptrm_check_for_SCAT.append(self.Data[s]['y_ptrm_check'][k])
          # If triangle is within the interval but started after the upper temperature bound, then one pTRM check is included
          # For example: if T_max=480, the traingle in 450 fall far, and it started at 500, then it is included
          # the ateration occured between 450 and 500, we dont know when.
          if  stop_scat_collect==False and \
             self.Data[s]['ptrm_checks_temperatures'][k] < pars["measurement_step_max"] and self.Data[s]['ptrm_checks_starting_temperatures'][k] > pars["measurement_step_max"] :
                x_ptrm_check_for_SCAT.append(self.Data[s]['x_ptrm_check'][k])
                y_ptrm_check_for_SCAT.append(self.Data[s]['y_ptrm_check'][k])
                stop_scat_collect=True
              
              
        # scat uses a different definistion":
        # use only pTRM that STARTED before the last temperatire step.
        
        x_ptrm_check_in_0_to_end=array(x_ptrm_check_in_0_to_end)  
        y_ptrm_check_in_0_to_end=array(y_ptrm_check_in_0_to_end)
        x_Arai_compare=array(x_Arai_compare)
        x_ptrm_check_in_start_to_end=array(x_ptrm_check_in_start_to_end)
        y_ptrm_check_in_start_to_end=array(y_ptrm_check_in_start_to_end)
        x_ptrm_check_for_SCAT=array(x_ptrm_check_for_SCAT)
        y_ptrm_check_for_SCAT=array(y_ptrm_check_for_SCAT)
                               
        DRATS=100*(abs(sum(x_ptrm_check_in_0_to_end-x_Arai_compare))/(x_Arai[end]))
        int_ptrm_n=len(x_ptrm_check_in_0_to_end)
        if int_ptrm_n > 0:
           pars['specimen_int_ptrm_n']=int_ptrm_n
           pars['specimen_drats']=DRATS
        else:
           pars['specimen_int_ptrm_n']=int_ptrm_n
           pars['specimen_drats']=-1

        #-------------------------------------------------
        # Tail check MD
        #-------------------------------------------------

        # collect tail check data"
        x_tail_check_start_to_end,y_tail_check_start_to_end=[],[]
        x_tail_check_for_SCAT,y_tail_check_for_SCAT=[],[]

        for k in range(len(self.Data[s]['tail_check_temperatures'])):
          if self.Data[s]['tail_check_temperatures'][k] in t_Arai:
              if self.Data[s]['tail_check_temperatures'][k]<=pars["measurement_step_max"] and self.Data[s]['tail_check_temperatures'][k] >=pars["measurement_step_min"]:
                   x_tail_check_start_to_end.append(self.Data[s]['x_tail_check'][k]) 
                   y_tail_check_start_to_end.append(self.Data[s]['y_tail_check'][k]) 
          if self.Data[s]['tail_check_temperatures'][k] >= pars["measurement_step_min"] and self.Data[s]['tail_checks_starting_temperatures'][k] <= pars["measurement_step_max"] :
                x_tail_check_for_SCAT.append(self.Data[s]['x_tail_check'][k])
                y_tail_check_for_SCAT.append(self.Data[s]['y_tail_check'][k])

                
        x_tail_check_start_to_end=array(x_tail_check_start_to_end)
        y_tail_check_start_to_end=array(y_tail_check_start_to_end)
        x_tail_check_for_SCAT=array(x_tail_check_for_SCAT)
        y_tail_check_for_SCAT=array(y_tail_check_for_SCAT)

        #-------------------------------------------------                     
        # Tail check : TO DO !
        pars['specimen_md']=-1  
        #-------------------------------------------------                     

        #-------------------------------------------------                     
        # Calculate the new 'beta box' parameter
        # all datapoints, pTRM checks, and tail-checks, should be inside a "beta box"
        # For definition of "beta box" see Shaar and Tauxe (2012)
        #-------------------------------------------------                     

        if self.accept_new_parameters['specimen_scat']==True or self.accept_new_parameters['specimen_scat'] in [1,"True","TRUE",'1']:
        
            pars["fail_arai_beta_box_scatter"]=False
            pars["fail_ptrm_beta_box_scatter"]=False
            pars["fail_tail_beta_box_scatter"]=False
            
            # best fit line 
            b=pars['specimen_b']
            cm_x=mean(array(x_Arai_segment))
            cm_y=mean(array(y_Arai_segment))
            a=cm_y-b*cm_x

            # lines with slope = slope +/- 2*(specimen_b_beta)

            if 'specimen_b_beta' not in self.accept_new_parameters.keys():
             print ("-E- ERROR: specimen_beta not in pmag_criteria file, cannot calculate 'beta box' scatter\n") 

            b_beta_threshold=self.accept_new_parameters['specimen_b_beta']

            two_sigma_beta_threshold=2*b_beta_threshold
            two_sigma_slope_threshold=abs(two_sigma_beta_threshold*b)
                 
            # a line with a  shallower  slope  (b + 2*beta*b) passing through the center of mass
            b1=b+two_sigma_slope_threshold
            a1=cm_y-b1*cm_x

            # bounding line with steeper  slope (b - 2*beta*b) passing through the center of mass
            b2=b-two_sigma_slope_threshold
            a2=cm_y-b2*cm_x

            # lower bounding line of the 'beta box'
            slop1=a1/((a2/b2))
            intercept1=a1

            # higher bounding line of the 'beta box'
            slop2=a2/((a1/b1))
            intercept2=a2       

            pars['specimen_scat_bounding_line_high']=[intercept2,slop2]
            pars['specimen_scat_bounding_line_low']=[intercept1,slop1]
            
            # check if the Arai data points are in the 'box'

            x_Arai_segment=array(x_Arai_segment)
            y_Arai_segment=array(y_Arai_segment)

            # the two bounding lines
            ymin=intercept1+x_Arai_segment*slop1
            ymax=intercept2+x_Arai_segment*slop2

            # arrays of "True" or "False"
            check_1=y_Arai_segment>ymax
            check_2=y_Arai_segment<ymin

            # check if at least one "True" 
            if (sum(check_1)+sum(check_2))>0:
             pars["fail_arai_beta_box_scatter"]=True
             #print "check, fail beta box"


            # check if the pTRM checks data points are in the 'box'

            # using x_ptrm_check_in_segment (defined above)
            # using y_ptrm_check_in_segment (defined above)


            if len(x_ptrm_check_for_SCAT) > 0:

              # the two bounding lines
              ymin=intercept1+x_ptrm_check_for_SCAT*slop1
              ymax=intercept2+x_ptrm_check_for_SCAT*slop2

              # arrays of "True" or "False"
              check_1=y_ptrm_check_for_SCAT>ymax
              check_2=y_ptrm_check_for_SCAT<ymin


              # check if at least one "True" 
              if (sum(check_1)+sum(check_2))>0:
                pars["fail_ptrm_beta_box_scatter"]=True
                #print "check, fail fail_ptrm_beta_box_scatter"
                
            # check if the tail checks data points are in the 'box'


            if len(x_tail_check_for_SCAT) > 0:

              # the two bounding lines
              ymin=intercept1+x_tail_check_for_SCAT*slop1
              ymax=intercept2+x_tail_check_for_SCAT*slop2

              # arrays of "True" or "False"
              check_1=y_tail_check_for_SCAT>ymax
              check_2=y_tail_check_for_SCAT<ymin


              # check if at least one "True" 
              if (sum(check_1)+sum(check_2))>0:
                pars["fail_tail_beta_box_scatter"]=True
                #print "check, fail fail_ptrm_beta_box_scatter"

            if pars["fail_tail_beta_box_scatter"] or pars["fail_ptrm_beta_box_scatter"] or pars["fail_arai_beta_box_scatter"]:
                  pars["specimen_scat"]="Fail"
            else:
                  pars["specimen_scat"]="Pass"
        else:
            pars["specimen_scat"]="N/A"
        #-------------------------------------------------  
        # Calculate the new FRAC parameter (Shaar and Tauxe, 2012).
        # also check that the 'gap' between consecutive measurements is less than 0.5(VDS)
        #
        #-------------------------------------------------  

        vector_diffs=self.Data[s]['vector_diffs']
        vector_diffs_segment=vector_diffs[zstart:zend]
        FRAC=sum(vector_diffs_segment)/self.Data[s]['vds']
        max_FRAC_gap=max(vector_diffs_segment/sum(vector_diffs_segment))

        pars['specimen_frac']=FRAC
        pars['specimen_gmax']=max_FRAC_gap

        #-------------------------------------------------  
        # Check if specimen pass Acceptance criteria
        #-------------------------------------------------  

        pars['specimen_fail_criteria']=[]
        for key in self.high_threshold_velue_list:
            if key in ['specimen_gmax','specimen_b_beta']:
                value=round(pars[key],2)
            elif key in ['specimen_dang','specimen_int_mad']:
                value=round(pars[key],1)
            else:
                value=pars[key]
                
            if value>float(self.accept_new_parameters[key]):
                pars['specimen_fail_criteria'].append(key)
        for key in self.low_threshold_velue_list:
            if key in ['specimen_f','specimen_fvds','specimen_frac','specimen_g','specimen_q']:
                value=round(pars[key],2)
            else: 
                value=pars[key]
            if value < float(self.accept_new_parameters[key]):
                pars['specimen_fail_criteria'].append(key)
        if 'specimen_scat' in pars.keys():
            if pars["specimen_scat"]=="Fail":
                pars['specimen_fail_criteria'].append('specimen_scat')
        if 'specimen_mad_scat' in pars.keys():
            if pars["specimen_mad_scat"]=="Fail":
                pars['specimen_fail_criteria'].append('specimen_mad_scat')

    
        #-------------------------------------------------                     
        # Calculate the direction of pTMRMS
        #-------------------------------------------------                     


        #-------------------------------------------------            
        # Calculate anistropy correction factor
        #-------------------------------------------------            

        if "AniSpec" in self.Data[s].keys():
           pars["AC_WARNING"]=""
           # if both aarm and atrm tensor axist, try first the aarm. if it fails use the atrm.
           if 'AARM' in self.Data[s]["AniSpec"].keys() and 'ATRM' in self.Data[s]["AniSpec"].keys():
               TYPES=['AARM','ATRM']
           else:
               TYPES=self.Data[s]["AniSpec"].keys()
           for TYPE in TYPES:
               red_flag=False
               S_matrix=zeros((3,3),'f')
               S_matrix[0,0]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s1']
               S_matrix[1,1]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s2']
               S_matrix[2,2]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s3']
               S_matrix[0,1]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s4']
               S_matrix[1,0]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s4']
               S_matrix[1,2]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s5']
               S_matrix[2,1]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s5']
               S_matrix[0,2]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s6']
               S_matrix[2,0]=self.Data[s]['AniSpec'][TYPE]['anisotropy_s6']

               #self.Data[s]['AniSpec']['anisotropy_type']=self.Data[s]['AniSpec']['anisotropy_type']
               self.Data[s]['AniSpec'][TYPE]['anisotropy_n']=int(float(self.Data[s]['AniSpec'][TYPE]['anisotropy_n']))

               this_specimen_f_type=self.Data[s]['AniSpec'][TYPE]['anisotropy_type']+"_"+"%i"%(int(self.Data[s]['AniSpec'][TYPE]['anisotropy_n']))
               
               Ftest_crit={} 
               Ftest_crit['ATRM_6']=  3.1059
               Ftest_crit['AARM_6']=  3.1059
               Ftest_crit['AARM_9']= 2.6848
               Ftest_crit['AARM_15']= 2.4558

               # threshold value for Ftest:
               
               if 'AniSpec' in self.Data[s].keys() and TYPE in self.Data[s]['AniSpec'].keys()\
                  and 'anisotropy_sigma' in  self.Data[s]['AniSpec'][TYPE].keys() \
                  and self.Data[s]['AniSpec'][TYPE]['anisotropy_sigma']!="":
                  # Calculate Ftest. If Ftest exceeds threshold value: set anistropy tensor to identity matrix
                   sigma=float(self.Data[s]['AniSpec'][TYPE]['anisotropy_sigma'])             
                   nf = 3*int(self.Data[s]['AniSpec'][TYPE]['anisotropy_n'])-6
                   F=calculate_ftest(S_matrix,sigma,nf)
                   #print s,"F",F
                   self.Data[s]['AniSpec'][TYPE]['ftest']=F
                   #print "s,sigma,nf,F,Ftest_crit[this_specimen_f_type]"
                   #print s,sigma,nf,F,Ftest_crit[this_specimen_f_type]
                   if self.accept_new_parameters['check_aniso_ftest']:
                       Ftest_threshold=Ftest_crit[this_specimen_f_type]
                       if self.Data[s]['AniSpec'][TYPE]['ftest'] < Ftest_crit[this_specimen_f_type]:
                           S_matrix=identity(3,'f')
                           pars["AC_WARNING"]=pars["AC_WARNING"]+"%s tensor fails F-test; "%(TYPE)
                           red_flag=True
                           
               else:
                   self.Data[s]['AniSpec'][TYPE]['anisotropy_sigma']=""
                   self.Data[s]['AniSpec'][TYPE]['ftest']=1e10
     
                
               if 'anisotropy_alt' in self.Data[s]['AniSpec'][TYPE].keys() and self.Data[s]['AniSpec'][TYPE]['anisotropy_alt']!="":
                   if float(self.Data[s]['AniSpec'][TYPE]['anisotropy_alt']) > float(self.accept_new_parameters['anisotropy_alt']):
                       S_matrix=identity(3,'f')
                       pars["AC_WARNING"]=pars["AC_WARNING"]+"%s tensor fails alteration check: %.1f%% > %.1f%%; "%(TYPE,float(self.Data[s]['AniSpec'][TYPE]['anisotropy_alt']),float(self.accept_new_parameters['anisotropy_alt']))
                       red_flag=True
               else:
                   self.Data[s]['AniSpec'][TYPE]['anisotropy_alt']=""

                   
               # if AARM passes, use the AARM.    
               if TYPE=='AARM' and red_flag==False:
                   break
               else:
                   pass
           TRM_anc_unit=array(pars['specimen_PCA_v1'])/sqrt(pars['specimen_PCA_v1'][0]**2+pars['specimen_PCA_v1'][1]**2+pars['specimen_PCA_v1'][2]**2)
           B_lab_unit=self.dir2cart([ self.Data[s]['Thellier_dc_field_phi'], self.Data[s]['Thellier_dc_field_theta'],1])
           #B_lab_unit=array([0,0,-1])
           Anisotropy_correction_factor=linalg.norm(dot(inv(S_matrix),TRM_anc_unit.transpose()))*norm(dot(S_matrix,B_lab_unit))
           pars["Anisotropy_correction_factor"]=Anisotropy_correction_factor
           pars["AC_specimen_int"]= pars["Anisotropy_correction_factor"] * float(pars["specimen_int"])
           
           pars["AC_anisotropy_type"]=self.Data[s]['AniSpec'][TYPE]["anisotropy_type"]
           pars["specimen_int_uT"]=float(pars["AC_specimen_int"])*1e6
           if TYPE=='AARM':
               if ":LP-AN-ARM" not in pars['magic_method_codes']:
                  pars['magic_method_codes']+=":LP-AN-ARM:AE-H:DA-AC-AARM"
                  pars['specimen_correction']='c'
                  pars['specimen_int_corr_anisotropy']=Anisotropy_correction_factor
           if TYPE=='ATRM':
               if ":LP-AN-TRM" not in pars['magic_method_codes']:
                  pars['magic_method_codes']+=":LP-AN-TRM:AE-H:DA-AC-ATRM"
                  pars['specimen_correction']='c' 
                  pars['specimen_int_corr_anisotropy']=Anisotropy_correction_factor

 
        else:
           pars["Anisotropy_correction_factor"]=1.0
           pars["specimen_int_uT"]=float(pars["specimen_int"])*1e6
           pars["AC_WARNING"]="No anistropy correction"
           pars['specimen_correction']='u' 

        pars["specimen_int_corr_anisotropy"]=pars["Anisotropy_correction_factor"]   
        #-------------------------------------------------                    
        # NLT and anisotropy correction together in one equation
        # See Shaar et al (2010), Equation (3)
        #-------------------------------------------------

        if 'NLT_parameters' in self.Data[s].keys():

           alpha=self.Data[s]['NLT_parameters']['tanh_parameters'][0][0]
           beta=self.Data[s]['NLT_parameters']['tanh_parameters'][0][1]
           b=float(pars["specimen_b"])
           Fa=pars["Anisotropy_correction_factor"]

           if ((abs(b)*Fa)/alpha) <1.0:
               Banc_NLT=math.atanh( ((abs(b)*Fa)/alpha) ) / beta
               pars["NLTC_specimen_int"]=Banc_NLT
               pars["specimen_int_uT"]=Banc_NLT*1e6

               if "AC_specimen_int" in pars.keys():
                   pars["NLT_specimen_correction_factor"]=Banc_NLT/float(pars["AC_specimen_int"])
               else:                       
                   pars["NLT_specimen_correction_factor"]=Banc_NLT/float(pars["specimen_int"])
               if ":LP-TRM" not in pars['magic_method_codes']:
                  pars['magic_method_codes']+=":LP-TRM:DA-NL"
               pars['specimen_correction']='c' 

           else:
               print ("-W- WARNING: problematic NLT mesurements for specimens %s. Cant do NLT calculation. check data\n"%s)
               pars["NLT_specimen_correction_factor"]=-1
        else:
           pars["NLT_specimen_correction_factor"]=-1

        #-------------------------------------------------                    
        # Calculate the final result with cooling rate correction
        #-------------------------------------------------

        pars["specimen_int_corr_cooling_rate"]=-999
        if 'cooling_rate_data' in self.Data[s].keys():
            if 'CR_correction_factor' in self.Data[s]['cooling_rate_data'].keys():
                if self.Data[s]['cooling_rate_data']['CR_correction_factor'] != -1 and self.Data[s]['cooling_rate_data']['CR_correction_factor'] !=-999:
                    pars["specimen_int_corr_cooling_rate"]=self.Data[s]['cooling_rate_data']['CR_correction_factor']
                    pars['specimen_correction']='c'
                    pars["specimen_int_uT"]=pars["specimen_int_uT"]*pars["specimen_int_corr_cooling_rate"]
                    if ":DA-CR" not in pars['magic_method_codes']:
                      pars['magic_method_codes']+=":DA-CR"
                    if 'CR_correction_factor_flag' in self.Data[s]['cooling_rate_data'].keys() \
                       and self.Data[s]['cooling_rate_data']['CR_correction_factor_flag']!="calculated":
                        pars["CR_WARNING"]="inferred cooling rate correction"
                    
                
        else:
            pars["CR_WARNING"]="no cooling rate correction"
            
##            sample=self.Data_hierarchy['specimens'][self.s]
##            if sample in Data_info["er_samples"]:
##                if 'sample_type' in Data_info["er_samples"][sample].keys():
##                    if Data_info["er_samples"][sample]['sample_type'] in ["Baked Clay","Baked Mud",

        print "pars before importing spd:  ", pars
        import spd
        print "imported spd"
        specimen = self.s
        pint_pars = spd.PintPars(self.Data,specimen,tmin,tmax)
        print "about to 'calculate all statistics'"
        print "pint_pars.calculate_all_statistics()"
        pint_pars.calculate_all_statistics()
        print "giraffes are awesome"
#        return(pars) original
        return pint_pars
    

    def get_default_criteria(self):
      #------------------------------------------------
      # read criteria file
      # Format is as pmag_criteria.txt
      #------------------------------------------------
      print "calling get_default_criteria()"

      self.criteria_list=['specimen_int_n','specimen_int_ptrm_n','specimen_f','specimen_fvds','specimen_frac','specimen_gmax','specimen_b_beta',
                     'specimen_dang','specimen_drats','specimen_int_mad','specimen_md','specimen_g','specimen_q']
      self.high_threshold_velue_list=['specimen_gmax','specimen_b_beta','specimen_dang','specimen_drats','specimen_int_mad','specimen_md']
      self.low_threshold_velue_list=['specimen_int_n','specimen_int_ptrm_n','specimen_f','specimen_fvds','specimen_frac','specimen_g','specimen_q']

      accept_new_parameters_null={}
      accept_new_parameters_default={}
      #  make a list of default parameters

      accept_new_parameters_default['specimen_int_n']=3
      accept_new_parameters_default['specimen_int_ptrm_n']=2
      accept_new_parameters_default['specimen_f']=0.
      accept_new_parameters_default['specimen_fvds']=0.
      accept_new_parameters_default['specimen_frac']=0.8
      accept_new_parameters_default['specimen_gmax']=0.6
      accept_new_parameters_default['specimen_b_beta']=0.1
      accept_new_parameters_default['specimen_dang']=100000
      accept_new_parameters_default['specimen_drats']=100000
      accept_new_parameters_default['specimen_int_mad']=5
      accept_new_parameters_default['specimen_md']=100000
      accept_new_parameters_default['specimen_g']=0
      accept_new_parameters_default['specimen_q']=0
      accept_new_parameters_default['specimen_scat']=True

      accept_new_parameters_default['sample_int_n']=3
      accept_new_parameters_default['sample_int_n_outlier_check']=6

      # anistropy criteria
      accept_new_parameters_default['anisotropy_alt']=10
      accept_new_parameters_default['check_aniso_ftest']=True


      # Sample mean calculation type 
      accept_new_parameters_default['sample_int_stdev_opt']=True
      accept_new_parameters_default['sample_int_bs']=False
      accept_new_parameters_default['sample_int_bs_par']=False

      # STDEV-OPT  
      accept_new_parameters_default['sample_int_sigma_uT']=6
      accept_new_parameters_default['sample_int_sigma_perc']=10
      accept_new_parameters_default['sample_aniso_threshold_perc']=1000000
      accept_new_parameters_default['sample_int_interval_uT']=10000
      accept_new_parameters_default['sample_int_interval_perc']=10000

      # BS  
      accept_new_parameters_default['sample_int_BS_68_uT']=10000
      accept_new_parameters_default['sample_int_BS_68_perc']=10000
      accept_new_parameters_default['sample_int_BS_95_uT']=10000
      accept_new_parameters_default['sample_int_BS_95_perc']=10000
      accept_new_parameters_default['specimen_int_max_slope_diff']=10000

      #    
      # NULL  
      for key in ( accept_new_parameters_default.keys()):
          accept_new_parameters_null[key]=accept_new_parameters_default[key]
      accept_new_parameters_null['sample_int_stdev_opt']=False
      accept_new_parameters_null['specimen_frac']=0
      accept_new_parameters_null['specimen_gmax']=10000
      accept_new_parameters_null['specimen_b_beta']=10000
      accept_new_parameters_null['specimen_int_mad']=100000
      accept_new_parameters_null['specimen_scat']=False
      accept_new_parameters_null['specimen_int_ptrm_n']=0
      accept_new_parameters_null['anisotropy_alt']=1e10
      accept_new_parameters_null['check_aniso_ftest']=True
      accept_new_parameters_default['sample_aniso_threshold_perc']=1000000

      accept_new_parameters_null['sample_int_sigma_uT']=0
      accept_new_parameters_null['sample_int_sigma_perc']=0
      accept_new_parameters_null['sample_int_n_outlier_check']=100000

      
      #print accept_new_parameters_default
        
      # A list of all acceptance criteria used by program
      accept_specimen_keys=['specimen_int_n','specimen_int_ptrm_n','specimen_f','specimen_fvds','specimen_frac','specimen_gmax','specimen_b_beta','specimen_dang','specimen_drats','specimen_int_mad','specimen_md']
      accept_sample_keys=['sample_int_n','sample_int_sigma_uT','sample_int_sigma_perc','sample_aniso_threshold_perc','sample_int_interval_uT','sample_int_interval_perc']
      
      #self.accept_new_parameters_null=accept_new_parameters_null
      return(accept_new_parameters_default,accept_new_parameters_null)
      #print accept_new_parameters_default
      #print "yes"
    
      
    def get_data(self):
      print "calling get_data()"
      print "self", self
      print "self.Data", self.Data
      print "magic file:", self.magic_file

      def tan_h(x, a, b):
          print "calling tan_h in get_data()"
          return a*tanh(b*x)
    

      #self.dir_pathes=self.WD


      #------------------------------------------------
      # Read magic measurement file and sort to blocks
      #------------------------------------------------

      # All data information is stored in Data[specimen]={}
      Data={}
      Data_hierarchy={}
      Data_hierarchy['samples']={}
      Data_hierarchy['specimens']={}

      # add dir to dir pathes for interpterer:
      if self.WD not in self.MagIC_directories_list:
          self.MagIC_directories_list.append(self.WD)
      #for dir_path in self.dir_pathes:
      #print "start Magic read %s " %self.magic_file
      try:
          meas_data,file_type=self.magic_read(self.magic_file)
      except:
          print "-E- ERROR: Cant read magic_measurement.txt file. File is corrupted."
          return {},{}

      #print "done Magic read %s " %self.magic_file

      print("-I- Read magic file  %s\n"%self.magic_file)

      # get list of unique specimen names
      
      CurrRec=[]
      #print "get sids"
      sids=self.get_specs(meas_data) # samples ID's
      #print "done get sids"

      #print "initialize blocks"
      
      for s in sids:

          if s not in Data.keys():
              Data[s]={}
              Data[s]['datablock']=[]
              Data[s]['trmblock']=[]
              Data[s]['zijdblock']=[]
          #zijdblock,units=pmag.find_dmag_rec(s,meas_data)
          #Data[s]['zijdblock']=zijdblock


      #print "done initialize blocks"

      #print "sorting meas data"
          
      for rec in meas_data:
          s=rec["er_specimen_name"]
          Data[s]['T_or_MW']="T"
          sample=rec["er_sample_name"]

          if  "LP-PI-M" in rec["magic_method_codes"]:
             Data[s]['T_or_MW']="MW"
          else:
             Data[s]['T_or_MW']="T"

          if "magic_method_codes" not in rec.keys():
              rec["magic_method_codes"]=""
          #methods=rec["magic_method_codes"].split(":")
          if "LP-PI-TRM" in rec["magic_method_codes"] or "LP-PI-M" in rec["magic_method_codes"]:
              Data[s]['datablock'].append(rec)
              # identify the lab DC field
              if ("LT-PTRM-I" in rec["magic_method_codes"] and 'LP-TRM' not in rec["magic_method_codes"] ) or "LT-PMRM-I" in rec["magic_method_codes"]:
                  Data[s]['Thellier_dc_field_uT']=float(rec["treatment_dc_field"])
                  Data[s]['Thellier_dc_field_phi']=float(rec['treatment_dc_field_phi'])
                  Data[s]['Thellier_dc_field_theta']=float(rec['treatment_dc_field_theta'])

                  
                
          if "LP-TRM" in rec["magic_method_codes"]:
              Data[s]['trmblock'].append(rec)

          if "LP-AN-TRM" in rec["magic_method_codes"]:
              if 'atrmblock' not in Data[s].keys():
                Data[s]['atrmblock']=[]
              Data[s]['atrmblock'].append(rec)


          if "LP-AN-ARM" in rec["magic_method_codes"]:
              if 'aarmblock' not in Data[s].keys():
                Data[s]['aarmblock']=[]
              Data[s]['aarmblock'].append(rec)

          if "LP-CR-TRM" in rec["magic_method_codes"]:
              if 'crblock' not in Data[s].keys():
                Data[s]['crblock']=[]
              Data[s]['crblock'].append(rec)

          #---- Zijderveld block

          EX=["LP-AN-ARM","LP-AN-TRM","LP-ARM-AFD","LP-ARM2-AFD","LP-TRM-AFD","LP-TRM","LP-TRM-TD","LP-X"] # list of excluded lab protocols
          #INC=["LT-NO","LT-AF-Z","LT-T-Z", "LT-M-Z", "LP-PI-TRM-IZ", "LP-PI-M-IZ"]
          INC=["LT-NO","LT-T-Z","LT-M-Z"]
          methods=rec["magic_method_codes"].split(":")
          for i in range (len(methods)):
               methods[i]=methods[i].strip()
          if 'measurement_flag' not in rec.keys(): rec['measurement_flag']='g'
          skip=1
          for meth in methods:
               if meth in INC:
                   skip=0
          for meth in EX:
               if meth in methods:skip=1
          if skip==0:
             if  'treatment_temp' in rec.keys():
                 tr = float(rec["treatment_temp"])
             elif "treatment_mw_power" in rec.keys():
                 tr = float(rec["treatment_mw_power"])
                 
             if "LP-PI-TRM-IZ" in methods or "LP-PI-M-IZ" in methods:  # looking for in-field first thellier or microwave data - otherwise, just ignore this
                 ZI=0
             else:
                 ZI=1
             Mkeys=['measurement_magnitude','measurement_magn_moment','measurement_magn_volume','measurement_magn_mass']
             if tr !="":
                 dec,inc,int = "","",""
                 if "measurement_dec" in rec.keys() and rec["measurement_dec"] != "":
                     dec=float(rec["measurement_dec"])
                 if "measurement_inc" in rec.keys() and rec["measurement_inc"] != "":
                     inc=float(rec["measurement_inc"])
                 for key in Mkeys:
                     if key in rec.keys() and rec[key]!="":int=float(rec[key])
                 if 'magic_instrument_codes' not in rec.keys():rec['magic_instrument_codes']=''
                 #datablock.append([tr,dec,inc,int,ZI,rec['measurement_flag'],rec['magic_instrument_codes']])
                 if Data[s]['T_or_MW']=="T":
                     if tr==0.: tr=273.
                 Data[s]['zijdblock'].append([tr,dec,inc,int,ZI,rec['measurement_flag'],rec['magic_instrument_codes']])
                 #print methods

       
          if sample not in Data_hierarchy['samples'].keys():
              Data_hierarchy['samples'][sample]=[]
          if s not in Data_hierarchy['samples'][sample]:
              Data_hierarchy['samples'][sample].append(s)

          Data_hierarchy['specimens'][s]=sample


          
      #print "done sorting meas data"
      
      self.specimens=Data.keys()
      self.s = self.specimens[0]  # LORI WEIRD ADDITION
      self.specimens.sort()

      
      #------------------------------------------------
      # Read anisotropy file from rmag_anisotropy.txt
      #------------------------------------------------

      #if self.WD != "":
      rmag_anis_data=[]
      results_anis_data=[]
      try:
          rmag_anis_data,file_type=self.magic_read(self.WD+'/rmag_anisotropy.txt')
          print( "-I- Anisotropy data read  %s/from rmag_anisotropy.txt\n"%self.WD)
      except:
          print("-W- WARNING cant find rmag_anisotropy in working directory\n")

      try:
          results_anis_data,file_type=self.magic_read(self.WD+'/rmag_results.txt')
          print( "-I- Anisotropy data read  %s/from rmag_anisotropy.txt\n"%self.WD)
          
      except:
          print("-W- WARNING cant find rmag_anisotropy in working directory\n")

          
      for AniSpec in rmag_anis_data:
          s=AniSpec['er_specimen_name']

          if s not in Data.keys():
              print("-W- WARNING: specimen %s in rmag_anisotropy.txt but not in magic_measurement.txt. Check it !\n"%s)
              continue
          if 'AniSpec' in Data[s].keys():
              print("-W- WARNING: more than one anisotropy data for specimen %s !\n"%s)
          TYPE=AniSpec['anisotropy_type']
          if 'AniSpec' not in Data[s].keys():
              Data[s]['AniSpec']={}
          Data[s]['AniSpec'][TYPE]=AniSpec
        
      for AniSpec in results_anis_data:
          s=AniSpec['er_specimen_names']
          if s not in Data.keys():
              print("-W- WARNING: specimen %s in rmag_results.txt but not in magic_measurement.txt. Check it !\n"%s)
              continue
          TYPE=AniSpec['anisotropy_type']         
          if 'AniSpec' in Data[s].keys() and TYPE in  Data[s]['AniSpec'].keys():
              Data[s]['AniSpec'][TYPE].update(AniSpec)
              if 'result_description' in AniSpec.keys():
                result_description=AniSpec['result_description'].split(";")
                for description in result_description:
                    if "Critical F" in description:
                       desc=description.split(":")
                       Data[s]['AniSpec'][TYPE]['anisotropy_F_crit']=float(desc[1])
                
                          
      #------------------------------------------------
      # Calculate Non Linear TRM parameters
      # Following Shaar et al. (2010):
      #
      # Procedure:
      #
      # A) If there are only 2 NLT measurement: C
      #
      #   Cant do NLT correctio procedure (few data points).
      #   Instead, check the different in the ratio (M/B) in the two measurements.
      #   slop_diff = max(first slope, second slope)/min(first slope, second slope)
      #   if: 1.1 > slop_diff > 1.05 : WARNING
      #   if: > slop_diff > 1s.1 : WARNING
      #
      # B) If there are at least 3 NLT measurement:
      #
      # 1) Read the NLT measurement file
      #   If there is no baseline measurement in the NLT experiment:
      #    then take the baseline from the zero-field step of the IZZI experiment.
      #
      # 2) Fit tanh function of the NLT measurement normalized by M[oven field]
      #   M/M[oven field] = alpha * tanh (beta*B)
      #   alpha and beta are used for the Banc calculation using equation (3) in Shaar et al. (2010):
      #   Banc= tanh^-1[(b*Fa)/alpha]/beta where Fa  is anistropy correction factor and 'b' is the Arai plot slope.
      #
      # 3) If best fit function algorithm does not converge, check NLT data using option (A) above.
      #    If 
      #
      #------------------------------------------------



      # Searching and sorting NLT Data 
      #print "searching NLT data"

      for s in self.specimens:
          datablock = Data[s]['datablock']
          trmblock = Data[s]['trmblock']

          if len(trmblock)<2:
              continue

          B_NLT,M_NLT=[],[]

          # find temperature of NLT acquisition
          NLT_temperature=float(trmblock[0]['treatment_temp'])
          
                 
          # search for Blab used in the IZZI experiment (need it for the following calculation)
          found_labfield=False  
          for rec in datablock:  
              if float(rec['treatment_dc_field'])!=0:
                  labfield=float(rec['treatment_dc_field'])
                  found_labfield=True
                  break
          if not found_labfield:
              continue

          # collect the data from trmblock
          M_baseline=0.
          for rec in trmblock:

              # if there is a baseline in TRM block, then use it 
              if float(rec['treatment_dc_field'])==0:
                  M_baseline=float(rec['measurement_magn_moment'])
              B_NLT.append(float(rec['treatment_dc_field']))
              M_NLT.append(float(rec['measurement_magn_moment']))

          # collect more data from araiblock


          for rec in datablock:
              if float(rec['treatment_temp'])==NLT_temperature and float(rec['treatment_dc_field']) !=0:
                  B_NLT.append(float(rec['treatment_dc_field']))
                  M_NLT.append(float(rec['measurement_magn_moment']))
                  
    
          # If cnat find baseline in trm block
          #  search for baseline in the Data block. 
          if M_baseline==0:
              m_tmp=[]
              for rec in datablock:
                  if float(rec['treatment_temp'])==NLT_temperature and float(rec['treatment_dc_field'])==0:
                     m_tmp.append(float(rec['measurement_magn_moment']))
                     print("-I- Found basleine for NLT measurements in datablock, specimen %s\n"%s)         
              if len(m_tmp)>0:
                  M_baseline=mean(m_tmp)
              

          ####  Ron dont delete it ### print "-I- Found %i NLT datapoints for specimen %s: B="%(len(B_NLT),s),array(B_NLT)*1e6
    
          #substitute baseline
          M_NLT=array(M_NLT)-M_baseline
          B_NLT=array(B_NLT)  
          # calculate M/B ratio for each step, and compare them
          # If cant do NLT correction: check a difference in M/B ratio
          # > 5% : WARNING
          # > 10%: ERROR           

          slopes=M_NLT/B_NLT

          if len(trmblock)==2:
              if max(slopes)/min(slopes)<1.05:
                  print("-I- 2 NLT measurement for specimen %s. [max(M/B)/ [min(M/B)] < 1.05.\n"%s)         
              elif max(slopes)/min(slopes)<1.1:
                  print("-W- WARNING: 2 NLT measurement for specimen %s. [max(M/B)]/ [min(M/B)] is %.2f  (   > 1.05 and  < 1.1 ). More NLT mrasurements may be required.\n" %(s,max(slopes)/min(slopes)))
                  #print("-I- NLT meaurements specime %s: B,M="%s,B_NLT,M_NLT)
              else:
                  print("-E- ERROR: 2 NLT measurement for specimen %s. [max(M/B)]/ [min(M/B)] is %.2f  ( > 1.1 ). More NLT mrasurements may be required  !\n" %(s,max(slopes)/min(slopes)))
                  #print("-I- NLT meaurements specime %s: B,M="%s,B_NLT,M_NLT)
                  
          # NLT procedure following Shaar et al (2010)        
          
          if len(trmblock)>2:
              B_NLT=append([0.],B_NLT)
              M_NLT=append([0.],M_NLT)
              
              try:
                  #print s,B_NLT, M_NLT    
                  # First try to fit tanh function (add point 0,0 in the begining)
                  alpha_0=max(M_NLT)
                  beta_0=2e4
                  popt, pcov = curve_fit(tan_h, B_NLT, M_NLT,p0=(alpha_0,beta_0))
                  M_lab=popt[0]*math.tanh(labfield*popt[1])

                  # Now  fit tanh function to the normalized curve
                  M_NLT_norm=M_NLT/M_lab
                  popt, pcov = curve_fit(tan_h, B_NLT, M_NLT_norm,p0=(popt[0]/M_lab,popt[1]))
                  Data[s]['NLT_parameters']={}
                  Data[s]['NLT_parameters']['tanh_parameters']=(popt, pcov)
                  Data[s]['NLT_parameters']['B_NLT']=B_NLT
                  Data[s]['NLT_parameters']['M_NLT_norm']=M_NLT_norm
                  
                  print("-I-  tanh parameters for specimen %s were calculated sucsessfuly\n"%s)
                                  
              except RuntimeError:
                  print( "-W- WARNING: Cant fit tanh function to NLT data specimen %s. Ignore NLT data for specimen %s. Instead check [max(M/B)]/ [min(M/B)] \n"%(s,s))
                  #print "-I- NLT meaurements specime %s: B,M="%s,B_NLT,M_NLT
                  
                  # Cant do NLT correction. Instead, check a difference in M/B ratio
                  # The maximum difference allowd is 5%
                  # if difference is larger than 5%: WARNING            
                  
                  if max(slopes)/min(slopes)<1.05:
                      print("-I- 2 NLT measurement for specimen %s. [max(M/B)/ [min(M/B)] < 1.05.\n"%s)         
                  elif max(slopes)/min(slopes)<1.1:
                      print("-W- WARNING: 2 NLT measurement for specimen %s. [max(M/B)]/ [min(M/B)] is %.2f  (   > 1.05 and  < 1.1 ). More NLT mrasurements may be required.\n" %(s,max(slopes)/min(slopes)))
                      #print "-I- NLT meaurements specime %s: B,M="%s,B_NLT,M_NLT
                  else:
                      print("-E- ERROR: 2 NLT measurement for specimen %s. [max(M/B)]/ [min(M/B)] is %.2f  ( > 1.1 ). More NLT mrasurements may be required  !\n" %(s,max(slopes)/min(slopes)))
                      #print "-I- NLT meaurements specime %s: B,M="%s,B_NLT,M_NLT
                  
      #print "done searching NLT data"
              
      print("-I- Done calculating non linear TRM parameters for all specimens\n")


      #------------------------------------------------
      # Calculate cooling rate experiments
      #
      #
      #
      #
      #
      #------------------------------------------------

      for s in self.specimens:
          datablock = Data[s]['datablock']
          trmblock = Data[s]['trmblock']
          if 'crblock' in Data[s].keys():
              if len(Data[s]['crblock'])<3:
                  del Data[s]['crblock']
                  continue

              sample=Data_hierarchy['specimens'][s]
              # in MagIC format that cooling rate is in K/My
##              try:
##                  ancient_cooling_rate=float(self.Data_info["er_samples"][sample]['sample_cooling_rate'])
##                  ancient_cooling_rate=ancient_cooling_rate/(1e6*365*24*60) # change to K/minute
##              except:
##                  print("-W- Cant find ancient cooling rate estimation for sample %s"%sample)
##                  continue                  
              try:
                  ancient_cooling_rate=float(self.Data_info["er_samples"][sample]['sample_cooling_rate'])
                  ancient_cooling_rate=ancient_cooling_rate/(1e6*365*24*60) # change to K/minute
              except:
                  print("-W- Cant find ancient cooling rate estimation for sample %s"%sample)
                  continue
              self.Data_info["er_samples"]
              cooling_rate_data={}
              cooling_rate_data['pairs']=[]
              cooling_rates_list=[]
              cooling_rate_data['alteration_check']=[]
              for rec in Data[s]['crblock']:
                  magic_method_codes=rec['magic_method_codes'].strip(' ').strip('\n').split(":")
                  measurement_description=rec['measurement_description'].strip(' ').strip('\n').split(":")
                  if "LT-T-Z" in magic_method_codes:
                      cooling_rate_data['baseline']=float(rec['measurement_magn_moment'])
                      continue
                
                  index=measurement_description.index("K/min")
                  cooling_rate=float(measurement_description[index-1])
                  cooling_rates_list.append(cooling_rate)
                  moment=float(rec['measurement_magn_moment'])
                  if "LT-T-I" in magic_method_codes:
                      cooling_rate_data['pairs'].append([cooling_rate,moment])
                  if "LT-PTRM-I" in magic_method_codes:
                      cooling_rate_data['alteration_check']=[cooling_rate,moment]
              lab_cooling_rate=max(cooling_rates_list) 
              cooling_rate_data['lab_cooling_rate']= lab_cooling_rate                  

              #lab_cooling_rate = self.Data[self.s]['cooling_rate_data']['lab_cooling_rate']
              moments=[]
              lab_fast_cr_moments=[]
              lan_cooling_rates=[]
              for pair in cooling_rate_data['pairs']:
                    lan_cooling_rates.append(math.log(cooling_rate_data['lab_cooling_rate']/pair[0]))
                    moments.append(pair[1])
                    if pair[0]==cooling_rate_data['lab_cooling_rate']:
                        lab_fast_cr_moments.append(pair[1])
              #print s, cooling_rate_data['alteration_check']
              lan_cooling_rates.append(math.log(cooling_rate_data['lab_cooling_rate']/cooling_rate_data['alteration_check'][0]))
              lab_fast_cr_moments.append(cooling_rate_data['alteration_check'][1])
              moments.append(cooling_rate_data['alteration_check'][1])        

              lab_fast_cr_moment=mean(lab_fast_cr_moments)
              moment_norm=array(moments)/lab_fast_cr_moment
              (a,b)=polyfit(lan_cooling_rates, moment_norm, 1)
              #ancient_cooling_rate=0.41
              x0=math.log(lab_cooling_rate/ancient_cooling_rate)
              y0=a*x0+b
              MAX=max(lab_fast_cr_moments)
              MIN=min(lab_fast_cr_moments)
                      
              alteration_check_perc=100*abs((MAX-MIN)/mean(MAX,MIN))
              #print s,alteration_check_perc
              #print "--"
              cooling_rate_data['ancient_cooling_rate']=ancient_cooling_rate
              cooling_rate_data['CR_correction_factor']=-999
              cooling_rate_data['lan_cooling_rates']=lan_cooling_rates
              cooling_rate_data['moment_norm']=moment_norm
              cooling_rate_data['polyfit']=[a,b]
              cooling_rate_data['CR_correction_factor_flag']=""
              if y0<=1:
                  cooling_rate_data['CR_correction_factor_flag']=cooling_rate_data['CR_correction_factor_flag']+"bad CR measurement data "
                  cooling_rate_data['CR_correction_factor']=-999

              if alteration_check_perc>5:
                  cooling_rate_data['CR_correction_factor_flag']=cooling_rate_data['CR_correction_factor_flag']+"alteration < 5% "
                  cooling_rate_data['CR_correction_factor']=-999
              if y0>1 and alteration_check_perc<=5:    
                  cooling_rate_data['CR_correction_factor_flag']="calculated"
                  cooling_rate_data['CR_correction_factor']=1/(y0)
                  
              Data[s]['cooling_rate_data']= cooling_rate_data     

              
               
      # go over all specimens. if there is a specimen with no cooling rate data
      # use the mean cooling rate corretion of the othr specimens from the same sample
      # this cooling rate correction is flagges as "inferred"

      for sample in Data_hierarchy['samples'].keys():
          CR_corrections=[]
          for s in Data_hierarchy['samples'][sample]:
              if 'cooling_rate_data' in Data[s].keys():
                  if 'CR_correction_factor' in Data[s]['cooling_rate_data'].keys():
                      if 'CR_correction_factor_flag' in Data[s]['cooling_rate_data'].keys():
                          if Data[s]['cooling_rate_data']['CR_correction_factor_flag']=='calculated':
                              CR_corrections.append(Data[s]['cooling_rate_data']['CR_correction_factor'])
          if len(CR_corrections) > 0:
              mean_CR_correction=mean(CR_corrections)
          else:
              mean_CR_correction=-1
          if mean_CR_correction != -1:
              for s in Data_hierarchy['samples'][sample]:
                  if 'cooling_rate_data' not in Data[s].keys():
                      Data[s]['cooling_rate_data']={}
                  if 'CR_correction_factor' not in Data[s]['cooling_rate_data'].keys() or\
                     Data[s]['cooling_rate_data']['CR_correction_factor_flag']!="calculated":
                        Data[s]['cooling_rate_data']['CR_correction_factor']=mean_CR_correction
                        Data[s]['cooling_rate_data']['CR_correction_factor_flag']="inferred"
              
      #------------------------------------------------
      # sort Arai block
      #------------------------------------------------

      #print "sort blocks to arai, zij. etc."

      for s in self.specimens:
        # collected the data
        datablock = Data[s]['datablock']
        zijdblock=Data[s]['zijdblock']

        if len(datablock) <4:
           print("-E- ERROR: skipping specimen %s, not enough measurements - moving forward \n"%s)
           del Data[s]
           sample=Data_hierarchy['specimens'][s]
           del Data_hierarchy['specimens'][s]
           Data_hierarchy['samples'][sample].remove(s)
           continue 
        araiblock,field=self.sortarai(datablock,s,0)


        Data[s]['araiblock']=araiblock
        Data[s]['pars']={}
        Data[s]['pars']['lab_dc_field']=field
        Data[s]['pars']['er_specimen_name']=s
        Data[s]['pars']['er_sample_name']=Data_hierarchy['specimens'][s]

        Data[s]['lab_dc_field']=field
        Data[s]['er_specimen_name']=s   
        Data[s]['er_sample_name']=Data_hierarchy['specimens'][s]
        
        first_Z=araiblock[0]
        #if len(first_Z)<3:
            #continue

        if len(araiblock[0])!= len(araiblock[1]):
           print( "-E- ERROR: unequal length of Z steps and I steps. Check specimen %s"% s)
           #continue

      # Fix zijderveld block for Thellier-Thellier protocol (II)
      # (take the vector subtruiction instead of the zerofield steps)
      #araiblock,field=self.sortarai(Data[s]['datablock'],s,0)
      #if "LP-PI-II" in Data[s]['datablock'][0]["magic_method_codes"] or "LP-PI-M-II" in Data[s]['datablock'][0]["magic_method_codes"] or "LP-PI-T-II" in Data[s]['datablock'][0]["magic_method_codes"]:
      #    for zerofield in araiblock[0]:
      #        Data[s]['zijdblock'].append([zerofield[0],zerofield[1],zerofield[2],zerofield[3],0,'g',""])
        if "LP-PI-II" in datablock[0]["magic_method_codes"] or "LP-PI-M-II" in datablock[0]["magic_method_codes"] or "LP-PI-T-II" in datablock[0]["magic_method_codes"]:
          for zerofield in araiblock[0]:
              Data[s]['zijdblock'].append([zerofield[0],zerofield[1],zerofield[2],zerofield[3],0,'g',""])


        #--------------------------------------------------------------
        # collect all zijderveld data to array and calculate VDS
        #--------------------------------------------------------------

        z_temperatures=[row[0] for row in zijdblock]
        zdata=[]
        vector_diffs=[]
        NRM=zijdblock[0][3]

        for k in range(len(zijdblock)):
            DIR=[zijdblock[k][1],zijdblock[k][2],zijdblock[k][3]/NRM]
            cart=self.dir2cart(DIR)
            zdata.append(array([cart[0],cart[1],cart[2]]))
            if k>0:
                vector_diffs.append(sqrt(sum((array(zdata[-2])-array(zdata[-1]))**2)))
        vector_diffs.append(sqrt(sum(array(zdata[-1])**2))) # last vector of the vds
        vds=sum(vector_diffs)  # vds calculation       
        zdata=array(zdata)
    
        Data[s]['vector_diffs']=array(vector_diffs)
        Data[s]['vds']=vds
        Data[s]['zdata']=zdata
        Data[s]['z_temp']=z_temperatures
        
      #--------------------------------------------------------------    
      # Rotate zijderveld plot
      #--------------------------------------------------------------

        DIR_rot=[]
        CART_rot=[]
        # rotate to be as NRM
        NRM_dir=self.cart2dir(Data[s]['zdata'][0])
         
        NRM_dec=NRM_dir[0]
        NRM_dir[0]=0
        CART_rot.append(self.dir2cart(NRM_dir))

        
        for i in range(1,len(Data[s]['zdata'])):
          DIR=self.cart2dir(Data[s]['zdata'][i])
          DIR[0]=DIR[0]-NRM_dec
          CART_rot.append(array(self.dir2cart(DIR)))
          #print array(dir2cart(DIR))
          
        CART_rot=array(CART_rot)
        Data[s]['zij_rotated']=CART_rot
        #--------------------------------------------------------------
        # collect all Arai plot data points to array 
        #--------------------------------------------------------------

        # collect Arai data points
        zerofields,infields=araiblock[0],araiblock[1]

        Data[s]['NRMS']=zerofields
        Data[s]['PTRMS']=infields
        
        x_Arai,y_Arai=[],[] # all the data points               
        t_Arai=[]
        steps_Arai=[]              

        #NRM=zerofields[0][3]
        infield_temperatures=[row[0] for row in infields]

        for k in range(len(zerofields)):                  
          index_infield=infield_temperatures.index(zerofields[k][0])
          x_Arai.append(infields[index_infield][3]/NRM)
          y_Arai.append(zerofields[k][3]/NRM)
          t_Arai.append(zerofields[k][0])
          if zerofields[k][4]==1:
            steps_Arai.append('ZI')
          else:
            steps_Arai.append('IZ')        
        x_Arai=array(x_Arai)
        y_Arai=array(y_Arai)
        #else:
        #    Data[s]['pars']['magic_method_codes']=""
        Data[s]['x_Arai']=x_Arai
        Data[s]['y_Arai']=y_Arai
        Data[s]['t_Arai']=t_Arai
        Data[s]['steps_Arai']=steps_Arai


        #--------------------------------------------------------------
        # collect all pTRM check to array 
        #--------------------------------------------------------------

        ptrm_checks = araiblock[2]
        zerofield_temperatures=[row[0] for row in zerofields]

        x_ptrm_check,y_ptrm_check,ptrm_checks_temperatures,=[],[],[]
        x_ptrm_check_starting_point,y_ptrm_check_starting_point,ptrm_checks_starting_temperatures=[],[],[]
        for k in range(len(ptrm_checks)):
          if ptrm_checks[k][0] in zerofield_temperatures:
            # find the starting point of the pTRM check:
            for i in range(len(datablock)):
                rec=datablock[i]                
                if "LT-PTRM-I" in rec['magic_method_codes'] and float(rec['treatment_temp'])==ptrm_checks[k][0]:
                    starting_temperature=(float(datablock[i-1]['treatment_temp']))

                    try:
                        index=t_Arai.index(starting_temperature)
                        x_ptrm_check_starting_point.append(x_Arai[index])
                        y_ptrm_check_starting_point.append(y_Arai[index])
                        ptrm_checks_starting_temperatures.append(starting_temperature)

                        index_zerofield=zerofield_temperatures.index(ptrm_checks[k][0])
                        x_ptrm_check.append(ptrm_checks[k][3]/NRM)
                        y_ptrm_check.append(zerofields[index_zerofield][3]/NRM)
                        ptrm_checks_temperatures.append(ptrm_checks[k][0])
                    except:
                        pass
                    
                # microwave
                if "LT-PMRM-I" in rec['magic_method_codes'] and float(rec['treatment_mw_power'])==ptrm_checks[k][0]:
                    starting_temperature=(float(datablock[i-1]['treatment_mw_power']))
                    
                    try:
                        index=t_Arai.index(starting_temperature)
                        x_ptrm_check_starting_point.append(x_Arai[index])
                        y_ptrm_check_starting_point.append(y_Arai[index])
                        ptrm_checks_starting_temperatures.append(starting_temperature)

                        index_zerofield=zerofield_temperatures.index(ptrm_checks[k][0])
                        x_ptrm_check.append(ptrm_checks[k][3]/NRM)
                        y_ptrm_check.append(zerofields[index_zerofield][3]/NRM)
                        ptrm_checks_temperatures.append(ptrm_checks[k][0])
                    except:
                        pass

                    
        x_ptrm_check=array(x_ptrm_check)  
        ptrm_check=array(y_ptrm_check)
        ptrm_checks_temperatures=array(ptrm_checks_temperatures)
        Data[s]['x_ptrm_check']=x_ptrm_check
        Data[s]['y_ptrm_check']=y_ptrm_check        
        Data[s]['ptrm_checks_temperatures']=ptrm_checks_temperatures
        Data[s]['x_ptrm_check_starting_point']=array(x_ptrm_check_starting_point)
        Data[s]['y_ptrm_check_starting_point']=array(y_ptrm_check_starting_point)               
        Data[s]['ptrm_checks_starting_temperatures']=array(ptrm_checks_starting_temperatures)
##        if len(ptrm_checks_starting_temperatures) != len(ptrm_checks_temperatures):
##            print s
##            print Data[s]['ptrm_checks_temperatures']
##            print Data[s]['ptrm_checks_starting_temperatures']
##            print "help"
            
        #--------------------------------------------------------------
        # collect tail checks 
        #--------------------------------------------------------------


        ptrm_tail = araiblock[3]
        #print ptrm_tail
        x_tail_check,y_tail_check,tail_check_temperatures=[],[],[]
        x_tail_check_starting_point,y_tail_check_starting_point,tail_checks_starting_temperatures=[],[],[]

        for k in range(len(ptrm_tail)):
          if ptrm_tail[k][0] in zerofield_temperatures:

            # find the starting point of the pTRM check:
            for i in range(len(datablock)):
                rec=datablock[i]                
                if "LT-PTRM-MD" in rec['magic_method_codes'] and float(rec['treatment_temp'])==ptrm_tail[k][0]:
                    starting_temperature=(float(datablock[i-1]['treatment_temp']))
                    try:

                        index=t_Arai.index(starting_temperature)
                        x_tail_check_starting_point.append(x_Arai[index])
                        y_tail_check_starting_point.append(y_Arai[index])
                        tail_checks_starting_temperatures.append(starting_temperature)

                        index_infield=infield_temperatures.index(ptrm_tail[k][0])
                        x_tail_check.append(infields[index_infield][3]/NRM)
                        y_tail_check.append(ptrm_tail[k][3]/NRM + zerofields[index_infield][3]/NRM)
                        tail_check_temperatures.append(ptrm_tail[k][0])

                        break
                    except:
                        pass


##              index_infield=infield_temperatures.index(ptrm_tail[k][0])
##              x_tail_check.append(infields[index_infield][3]/NRM)
##              y_tail_check.append(ptrm_tail[k][3]/NRM + zerofields[index_infield][3]/NRM)
##              tail_check_temperatures.append(ptrm_tail[k][0])

        x_tail_check=array(x_tail_check)  
        y_tail_check=array(y_tail_check)
        tail_check_temperatures=array(tail_check_temperatures)
        x_tail_check_starting_point=array(x_tail_check_starting_point)
        y_tail_check_starting_point=array(y_tail_check_starting_point)
        tail_checks_starting_temperatures=array(tail_checks_starting_temperatures)
        
        Data[s]['x_tail_check']=x_tail_check
        Data[s]['y_tail_check']=y_tail_check
        Data[s]['tail_check_temperatures']=tail_check_temperatures
        Data[s]['x_tail_check_starting_point']=x_tail_check_starting_point
        Data[s]['y_tail_check_starting_point']=y_tail_check_starting_point
        Data[s]['tail_checks_starting_temperatures']=tail_checks_starting_temperatures
        
##        #--------------------------------------------------------------
##        # collect tail checks 
##        #--------------------------------------------------------------
##
##
##        ptrm_tail = araiblock[3]
##        #print ptrm_tail
##        x_tail_check,y_tail_check=[],[]
##
##        for k in range(len(ptrm_tail)):                  
##          index_infield=infield_temperatures.index(ptrm_tail[k][0])
##          x_tail_check.append(infields[index_infield][3]/NRM)
##          y_tail_check.append(ptrm_tail[k][3]/NRM + zerofields[index_infield][3]/NRM)
##          
##
##        x_tail_check=array(x_tail_check)  
##        y_tail_check=array(y_tail_check)
##
##        Data[s]['x_tail_check']=x_tail_check
##        Data[s]['y_tail_check']=y_tail_check

      print("-I- number of specimens in this project directory: %i\n"%len(self.specimens))
      print("-I- number of samples in this project directory: %i\n"%len(Data_hierarchy['samples'].keys()))

      #print "done sort blocks to arai, zij. etc."
      print "returning Data, data_hierarchy.  This is the completion of self.get_data().  printing Data['0238x5721062']"
      print str(Data["0238x5721062"])[:500] + "...."
      print "done with get_data"
      return(Data,Data_hierarchy)


      
 # zebra.  end of get_data()

    #--------------------------------------------------------------    
    # Read all information file (er_locations, er_samples, er_sites, er_ages)
    #--------------------------------------------------------------
    def get_data_info(self):
        print "calling get_data_info()"
        Data_info={}
        data_er_samples={}
        data_er_ages={}
        data_er_sites={}

    
        # samples
        # read_magic_file takes 2 args
        # other_read_magic_file takes 4 args (including self) 
        def read_magic_file(path,sort_by_this_name):
            # called for er_ages, er_sites, er_samples
            print "Calling read_magic_file() in get_data_info"
            print path
            DATA={}
            fin=open(path,'rU')
            fin.readline()
            line=fin.readline()
            header=line.strip('\n').split('\t')
            for line in fin.readlines():
                tmp_data={}
                tmp_line=line.strip('\n').split('\t')
                for i in range(len(tmp_line)):
                    tmp_data[header[i]]=tmp_line[i]
                DATA[tmp_data[sort_by_this_name]]=tmp_data
            fin.close()        
            print "Data from read_magic_file in get_data info:  ", DATA
            return(DATA)
        
        try:
            data_er_samples=read_magic_file(self.WD+"/er_samples.txt",'er_sample_name')
        except:
            print "-W- Cant find er_sample.txt in project directory\n"
    
        try:
            data_er_sites=read_magic_file(self.WD+"/er_sites.txt",'er_site_name')
        except:
            print ("-W- Cant find er_sites.txt in project directory\n")

        try:
            data_er_ages=read_magic_file(self.WD+"/er_ages.txt",'er_sample_name')
        except:
            try:
                data_er_ages=read_magic_file(self.WD+"/er_ages.txt",'er_site_name')
            except:    
                print ("-W- Cant find er_ages in project directory\n")
    

        Data_info["er_samples"]=data_er_samples
        Data_info["er_sites"]=data_er_sites
        Data_info["er_ages"]=data_er_ages
        
        print "data_info"
        print str(Data_info)[:500]
        return(Data_info)

    #--------------------------------------------------------------    
    # Read previose interpretation from pmag_specimens.txt (if exist)
    #--------------------------------------------------------------
    
    def get_previous_interpretation(self):
        print "calling get_previous_interpretation()"
        try:
            print ("-I- Read pmag_specimens.txt for previouse interpretation")
            print "about to call self.read_magic_file()"
            prev_pmag_specimen=self.classy_read_magic_file(self.WD+"/pmag_specimens.txt",1,'er_specimen_name')
            #f
            print "successfully read pmag_specimens"
            # first delete all previous interpretation
            for sp in self.Data.keys():
                del self.Data[sp]['pars']
                self.Data[sp]['pars']={}
                self.Data[sp]['pars']['lab_dc_field']=self.Data[sp]['lab_dc_field']
                self.Data[sp]['pars']['er_specimen_name']=self.Data[sp]['er_specimen_name']   
                self.Data[sp]['pars']['er_sample_name']=self.Data[sp]['er_sample_name']

            self.Data_samples={}
        
            for specimen in prev_pmag_specimen.keys():
              tmin_kelvin=float(prev_pmag_specimen[specimen]['measurement_step_min'])
              tmax_kelvin=float(prev_pmag_specimen[specimen]['measurement_step_max'])
              if specimen not in self.redo_specimens.keys():
                self.redo_specimens[specimen]={}
              self.redo_specimens[specimen]['t_min']=float(tmin_kelvin)
              self.redo_specimens[specimen]['t_max']=float(tmax_kelvin)
              if specimen in self.Data.keys():
                  if tmin_kelvin not in self.Data[specimen]['t_Arai'] or tmax_kelvin not in self.Data[specimen]['t_Arai'] :
                      print ("-W- WARNING: cant fit temperature bounds in the redo file to the actual measurement. specimen %s\n"%specimen)
                  else:
                      try:
                          self.Data[specimen]['pars']=self.get_PI_parameters(specimen,float(tmin_kelvin),float(tmax_kelvin))
                          self.Data[specimen]['pars']['saved']=True
                          # write intrepretation into sample data
                          sample=self.Data_hierarchy['specimens'][specimen]
                          if sample not in self.Data_samples.keys():
                              self.Data_samples[sample]={}
                          self.Data_samples[sample][specimen]=self.Data[specimen]['pars']['specimen_int_uT']
                      except:
                          print ("-E- ERROR. Cant calculate PI paremeters for specimen %s using redo file. Check!"%(specimen))
              else:
                  print ("-W- WARNING: Cant find specimen %s from redo file in measurement file!\n"%specimen)
        
            try:
                self.s
            except:
                self.s=self.specimens[0]
                    
            self.pars=self.Data[self.s]['pars']
            self.clear_boxes()
            self.draw_figure(self.s)
            self.update_GUI_with_new_interpretation()
        except Exception as ex:
            print "exception: ", ex
            return
                    


#===========================================================
#  definitions inherited from pmag.py
#===========================================================
    
                
    def cart2dir(self,cart):
        """
        converts a direction to cartesian coordinates
        """
#        print "calling cart2dir(), not in anything"
        cart=array(cart)
        rad=pi/180. # constant to convert degrees to radians
        if len(cart.shape)>1:
            Xs,Ys,Zs=cart[:,0],cart[:,1],cart[:,2]
        else: #single vector
            Xs,Ys,Zs=cart[0],cart[1],cart[2]
        Rs=sqrt(Xs**2+Ys**2+Zs**2) # calculate resultant vector length
        Decs=(arctan2(Ys,Xs)/rad)%360. # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
        try:
            Incs=arcsin(Zs/Rs)/rad # calculate inclination (converting to degrees) # 
        except:
            print 'trouble in cart2dir' # most likely division by zero somewhere
            return zeros(3)
            
        return array([Decs,Incs,Rs]).transpose() # return the directions list


    def dir2cart(self,d):
#        print "calling dir2cart(), not in anything"
       # converts list or array of vector directions, in degrees, to array of cartesian coordinates, in x,y,z
        ints=ones(len(d)).transpose() # get an array of ones to plug into dec,inc pairs
        d=array(d)
        rad=pi/180.
        if len(d.shape)>1: # array of vectors
            decs,incs=d[:,0]*rad,d[:,1]*rad
            if d.shape[1]==3: ints=d[:,2] # take the given lengths
        else: # single vector
            decs,incs=array(d[0])*rad,array(d[1])*rad
            if len(d)==3: 
                ints=array(d[2])
            else:
                ints=array([1.])
        cart= array([ints*cos(decs)*cos(incs),ints*sin(decs)*cos(incs),ints*sin(incs)]).transpose()
        return cart


    def b_vdm(self,B,lat):
        """ 
        Converts field values in tesla to v(a)dm in Am^2
        """
        print "calling b_vdm()"
        not_called = """
        rad=pi/180.
        fact=((6.371e6)**3)*1e7 # changed radius of the earth from 3.367e6 3/12/2010
        colat=(90.-lat) * rad
        return fact*B/(sqrt(1+3*(cos(colat)**2)))"""

    def dohext(self,nf,sigma,s):
        """
        calculates hext parameters for nf, sigma and s
        """
        #
        print "calling dohext()"
        not_called="""
        if nf==-1:return hextpars 
        f=sqrt(2.*self.fcalc(2,nf))
        t2sum=0
        tau,Vdir=self.doseigs(s)
        for i in range(3): t2sum+=tau[i]**2
        chibar=(s[0]+s[1]+s[2])/3.
        hpars={}
        hpars['F_crit']='%s'%(self.fcalc(5,nf))
        hpars['F12_crit']='%s'%(self.fcalc(2,nf))
        hpars["F"]=0.4*(t2sum-3*chibar**2)/(sigma**2)
        hpars["F12"]=0.5*((tau[0]-tau[1])/sigma)**2
        hpars["F23"]=0.5*((tau[1]-tau[2])/sigma)**2
        hpars["v1_dec"]=Vdir[0][0]
        hpars["v1_inc"]=Vdir[0][1]
        hpars["v2_dec"]=Vdir[1][0]
        hpars["v2_inc"]=Vdir[1][1]
        hpars["v3_dec"]=Vdir[2][0]
        hpars["v3_inc"]=Vdir[2][1]
        hpars["t1"]=tau[0]
        hpars["t2"]=tau[1]
        hpars["t3"]=tau[2]
        hpars["e12"]=arctan((f*sigma)/(2*abs(tau[0]-tau[1])))*180./pi
        hpars["e23"]=arctan((f*sigma)/(2*abs(tau[1]-tau[2])))*180./pi
        hpars["e13"]=arctan((f*sigma)/(2*abs(tau[0]-tau[2])))*180./pi
        return hpars"""

    def doseigs(self,s):
        """
        convert s format for eigenvalues and eigenvectors
        """
    #
        print "calling doseigs()"
        not_called="""
        A=self.s2a(s) # convert s to a (see Tauxe 1998)
        tau,V=self.tauV(A) # convert to eigenvalues (t), eigenvectors (V)
        Vdirs=[]
        for v in V: # convert from cartesian to direction
            Vdir= self.cart2dir(v)
            if Vdir[1]<0:
                Vdir[1]=-Vdir[1]
                Vdir[0]=(Vdir[0]+180.)%360.
            Vdirs.append([Vdir[0],Vdir[1]])
        return tau,Vdirs"""


    def tauV(self,T):
        """
        gets the eigenvalues (tau) and eigenvectors (V) from matrix T
        """
        print "calling tauV()"
        not_called="""
        t,V,tr=[],[],0.
        ind1,ind2,ind3=0,1,2
        evalues,evectmps=linalg.eig(T)
        evectors=transpose(evectmps)  # to make compatible with Numeric convention
        for tau in evalues:
            tr+=tau
        if tr!=0:
            for i in range(3):
                evalues[i]=evalues[i]/tr
        else:
            return t,V
    # sort evalues,evectors
        t1,t2,t3=0.,0.,1.
        for k in range(3):
            if evalues[k] > t1: 
                t1,ind1=evalues[k],k 
            if evalues[k] < t3: 
                t3,ind3=evalues[k],k 
        for k in range(3):
            if evalues[k] != t1 and evalues[k] != t3: 
                t2,ind2=evalues[k],k
        V.append(evectors[ind1])
        V.append(evectors[ind2])
        V.append(evectors[ind3])
        t.append(t1)
        t.append(t2)
        t.append(t3)
        return t,V"""
    

    def s2a(self,s):
        """
         convert 6 element "s" list to 3,3 a matrix (see Tauxe 1998)
        """
        print "calling s2a()"
        not_called="""
        a=zeros((3,3,),'f') # make the a matrix
        for i in range(3):
            a[i][i]=s[i]
        a[0][1],a[1][0]=s[3],s[3]
        a[1][2],a[2][1]=s[4],s[4]
        a[0][2],a[2][0]=s[5],s[5]
        return a"""

    
    def fcalc(self,col,row):
        """
      looks up f from ftables F(row,col), where row is number of degrees of freedom - this is 95% confidence (p=0.05)
        """
    #
        print "calling fcalc()"
        not_called = """
        if row>200:row=200
        if col>20:col=20
        ftest=array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
    [1, 161.469, 199.493, 215.737, 224.5, 230.066, 234.001, 236.772, 238.949, 240.496, 241.838, 242.968, 243.88, 244.798, 245.26, 245.956, 246.422, 246.89, 247.36, 247.596, 248.068],
    [2, 18.5128, 18.9995, 19.1642, 19.2467, 19.2969, 19.3299, 19.3536, 19.371, 19.3852, 19.3963, 19.4043, 19.4122, 19.4186, 19.425, 19.4297, 19.4329, 19.4377, 19.4409, 19.4425, 19.4457],
    [3, 10.1278, 9.5522, 9.2767, 9.1173, 9.0133, 8.9408, 8.8868, 8.8452, 8.8124, 8.7857, 8.7635, 8.7446, 8.7287, 8.715, 8.7028, 8.6923, 8.683, 8.6745, 8.667, 8.6602],
    [4, 7.7087, 6.9444, 6.5915, 6.3882, 6.2561, 6.1631, 6.0943, 6.0411, 5.9988, 5.9644, 5.9359, 5.9117, 5.8912, 5.8733, 5.8578, 5.844, 5.8319, 5.8211, 5.8113, 5.8025],
    [5, 6.608, 5.7861, 5.4095, 5.1922, 5.0503, 4.9503, 4.8759, 4.8184, 4.7725, 4.735, 4.7039, 4.6777, 4.6552, 4.6358, 4.6187, 4.6038, 4.5904, 4.5785, 4.5679, 4.5581],
    [6, 5.9874, 5.1433, 4.757, 4.5337, 4.3874, 4.2838, 4.2067, 4.1468, 4.099, 4.06, 4.0275, 3.9999, 3.9764, 3.956, 3.9381, 3.9223, 3.9083, 3.8957, 3.8844, 3.8742],
    [7, 5.5914, 4.7374, 4.3469, 4.1204, 3.9715, 3.866, 3.787, 3.7257, 3.6767, 3.6366, 3.603, 3.5747, 3.5504, 3.5292, 3.5107, 3.4944, 3.4799, 3.4669, 3.4552, 3.4445],
    [8, 5.3177, 4.459, 4.0662, 3.8378, 3.6875, 3.5806, 3.5004, 3.4381, 3.3881, 3.3472, 3.313, 3.2839, 3.259, 3.2374, 3.2184, 3.2017, 3.1867, 3.1733, 3.1613, 3.1503],
    [9, 5.1174, 4.2565, 3.8626, 3.6331, 3.4817, 3.3738, 3.2928, 3.2296, 3.1789, 3.1373, 3.1025, 3.0729, 3.0475, 3.0255, 3.0061, 2.989, 2.9737, 2.96, 2.9476, 2.9365],
    [10, 4.9647, 4.1028, 3.7083, 3.4781, 3.3258, 3.2171, 3.1355, 3.0717, 3.0204, 2.9782, 2.9429, 2.913, 2.8872, 2.8648, 2.845, 2.8276, 2.812, 2.7981, 2.7855, 2.774],
    [11, 4.8443, 3.9823, 3.5875, 3.3567, 3.2039, 3.0946, 3.0123, 2.948, 2.8962, 2.8536, 2.8179, 2.7876, 2.7614, 2.7386, 2.7186, 2.7009, 2.6851, 2.6709, 2.6581, 2.6464],
    [12, 4.7472, 3.8853, 3.4903, 3.2592, 3.1059, 2.9961, 2.9134, 2.8486, 2.7964, 2.7534, 2.7173, 2.6866, 2.6602, 2.6371, 2.6169, 2.5989, 2.5828, 2.5684, 2.5554, 2.5436],
    [13, 4.6672, 3.8055, 3.4106, 3.1791, 3.0255, 2.9153, 2.8321, 2.7669, 2.7144, 2.6711, 2.6347, 2.6037, 2.5769, 2.5536, 2.5331, 2.5149, 2.4987, 2.4841, 2.4709, 2.4589],
    [14, 4.6001, 3.7389, 3.3439, 3.1122, 2.9582, 2.8477, 2.7642, 2.6987, 2.6458, 2.6021, 2.5655, 2.5343, 2.5073, 2.4837, 2.463, 2.4446, 2.4282, 2.4134, 2.4, 2.3879],
    [15, 4.543, 3.6824, 3.2874, 3.0555, 2.9013, 2.7905, 2.7066, 2.6408, 2.5877, 2.5437, 2.5068, 2.4753, 2.4481, 2.4244, 2.4034, 2.3849, 2.3683, 2.3533, 2.3398, 2.3275],
    [16, 4.494, 3.6337, 3.2389, 3.0069, 2.8524, 2.7413, 2.6572, 2.5911, 2.5377, 2.4935, 2.4564, 2.4247, 2.3973, 2.3733, 2.3522, 2.3335, 2.3167, 2.3016, 2.288, 2.2756],
    [17, 4.4513, 3.5916, 3.1968, 2.9647, 2.81, 2.6987, 2.6143, 2.548, 2.4943, 2.4499, 2.4126, 2.3807, 2.3531, 2.329, 2.3077, 2.2888, 2.2719, 2.2567, 2.2429, 2.2303],
    [18, 4.4139, 3.5546, 3.1599, 2.9278, 2.7729, 2.6613, 2.5767, 2.5102, 2.4563, 2.4117, 2.3742, 2.3421, 2.3143, 2.29, 2.2686, 2.2496, 2.2325, 2.2172, 2.2033, 2.1906],
    [19, 4.3808, 3.5219, 3.1274, 2.8951, 2.7401, 2.6283, 2.5435, 2.4768, 2.4227, 2.378, 2.3402, 2.308, 2.28, 2.2556, 2.2341, 2.2149, 2.1977, 2.1823, 2.1683, 2.1555],
    [20, 4.3512, 3.4928, 3.0984, 2.8661, 2.7109, 2.599, 2.514, 2.4471, 2.3928, 2.3479, 2.31, 2.2776, 2.2495, 2.2249, 2.2033, 2.184, 2.1667, 2.1511, 2.137, 2.1242],
    [21, 4.3248, 3.4668, 3.0725, 2.8401, 2.6848, 2.5727, 2.4876, 2.4205, 2.3661, 2.3209, 2.2829, 2.2504, 2.2222, 2.1975, 2.1757, 2.1563, 2.1389, 2.1232, 2.109, 2.096],
    [22, 4.3009, 3.4434, 3.0492, 2.8167, 2.6613, 2.5491, 2.4638, 2.3965, 2.3419, 2.2967, 2.2585, 2.2258, 2.1975, 2.1727, 2.1508, 2.1313, 2.1138, 2.098, 2.0837, 2.0707],
    [23, 4.2794, 3.4221, 3.028, 2.7955, 2.64, 2.5276, 2.4422, 2.3748, 2.3201, 2.2747, 2.2364, 2.2036, 2.1752, 2.1503, 2.1282, 2.1086, 2.091, 2.0751, 2.0608, 2.0476],
    [24, 4.2597, 3.4029, 3.0088, 2.7763, 2.6206, 2.5082, 2.4226, 2.3551, 2.3003, 2.2547, 2.2163, 2.1834, 2.1548, 2.1298, 2.1077, 2.088, 2.0703, 2.0543, 2.0399, 2.0267],
    [25, 4.2417, 3.3852, 2.9913, 2.7587, 2.603, 2.4904, 2.4047, 2.3371, 2.2821, 2.2365, 2.1979, 2.1649, 2.1362, 2.1111, 2.0889, 2.0691, 2.0513, 2.0353, 2.0207, 2.0075],
    [26, 4.2252, 3.369, 2.9752, 2.7426, 2.5868, 2.4741, 2.3883, 2.3205, 2.2655, 2.2197, 2.1811, 2.1479, 2.1192, 2.094, 2.0716, 2.0518, 2.0339, 2.0178, 2.0032, 1.9898],
    [27, 4.21, 3.3542, 2.9603, 2.7277, 2.5719, 2.4591, 2.3732, 2.3053, 2.2501, 2.2043, 2.1656, 2.1323, 2.1035, 2.0782, 2.0558, 2.0358, 2.0179, 2.0017, 1.987, 1.9736],
    [28, 4.196, 3.3404, 2.9467, 2.7141, 2.5581, 2.4453, 2.3592, 2.2913, 2.236, 2.1901, 2.1512, 2.1179, 2.0889, 2.0636, 2.0411, 2.021, 2.0031, 1.9868, 1.972, 1.9586],
    [29, 4.1829, 3.3276, 2.9341, 2.7014, 2.5454, 2.4324, 2.3463, 2.2783, 2.2229, 2.1768, 2.1379, 2.1045, 2.0755, 2.05, 2.0275, 2.0074, 1.9893, 1.973, 1.9582, 1.9446],
    [30, 4.1709, 3.3158, 2.9223, 2.6896, 2.5335, 2.4205, 2.3343, 2.2662, 2.2107, 2.1646, 2.1255, 2.0921, 2.0629, 2.0374, 2.0148, 1.9946, 1.9765, 1.9601, 1.9452, 1.9317],
    [31, 4.1597, 3.3048, 2.9113, 2.6787, 2.5225, 2.4094, 2.3232, 2.2549, 2.1994, 2.1531, 2.1141, 2.0805, 2.0513, 2.0257, 2.003, 1.9828, 1.9646, 1.9481, 1.9332, 1.9196],
    [32, 4.1491, 3.2945, 2.9011, 2.6684, 2.5123, 2.3991, 2.3127, 2.2444, 2.1888, 2.1425, 2.1033, 2.0697, 2.0404, 2.0147, 1.992, 1.9717, 1.9534, 1.9369, 1.9219, 1.9083],
    [33, 4.1392, 3.2849, 2.8915, 2.6589, 2.5027, 2.3894, 2.303, 2.2346, 2.1789, 2.1325, 2.0933, 2.0596, 2.0302, 2.0045, 1.9817, 1.9613, 1.943, 1.9264, 1.9114, 1.8977],
    [34, 4.13, 3.2759, 2.8826, 2.6499, 2.4936, 2.3803, 2.2938, 2.2253, 2.1696, 2.1231, 2.0838, 2.05, 2.0207, 1.9949, 1.972, 1.9516, 1.9332, 1.9166, 1.9015, 1.8877],
    [35, 4.1214, 3.2674, 2.8742, 2.6415, 2.4851, 2.3718, 2.2852, 2.2167, 2.1608, 2.1143, 2.0749, 2.0411, 2.0117, 1.9858, 1.9629, 1.9424, 1.924, 1.9073, 1.8922, 1.8784],
    [36, 4.1132, 3.2594, 2.8663, 2.6335, 2.4771, 2.3637, 2.2771, 2.2085, 2.1526, 2.1061, 2.0666, 2.0327, 2.0032, 1.9773, 1.9543, 1.9338, 1.9153, 1.8986, 1.8834, 1.8696],
    [37, 4.1055, 3.2519, 2.8588, 2.6261, 2.4696, 2.3562, 2.2695, 2.2008, 2.1449, 2.0982, 2.0587, 2.0248, 1.9952, 1.9692, 1.9462, 1.9256, 1.9071, 1.8904, 1.8752, 1.8613],
    [38, 4.0981, 3.2448, 2.8517, 2.619, 2.4625, 2.349, 2.2623, 2.1935, 2.1375, 2.0909, 2.0513, 2.0173, 1.9877, 1.9617, 1.9386, 1.9179, 1.8994, 1.8826, 1.8673, 1.8534],
    [39, 4.0913, 3.2381, 2.8451, 2.6123, 2.4558, 2.3422, 2.2555, 2.1867, 2.1306, 2.0839, 2.0442, 2.0102, 1.9805, 1.9545, 1.9313, 1.9107, 1.8921, 1.8752, 1.8599, 1.8459],
    [40, 4.0848, 3.2317, 2.8388, 2.606, 2.4495, 2.3359, 2.249, 2.1802, 2.124, 2.0773, 2.0376, 2.0035, 1.9738, 1.9476, 1.9245, 1.9038, 1.8851, 1.8682, 1.8529, 1.8389],
    [41, 4.0786, 3.2257, 2.8328, 2.6, 2.4434, 2.3298, 2.2429, 2.174, 2.1178, 2.071, 2.0312, 1.9971, 1.9673, 1.9412, 1.9179, 1.8972, 1.8785, 1.8616, 1.8462, 1.8321],
    [42, 4.0727, 3.2199, 2.8271, 2.5943, 2.4377, 2.324, 2.2371, 2.1681, 2.1119, 2.065, 2.0252, 1.991, 1.9612, 1.935, 1.9118, 1.8909, 1.8722, 1.8553, 1.8399, 1.8258],
    [43, 4.067, 3.2145, 2.8216, 2.5888, 2.4322, 2.3185, 2.2315, 2.1625, 2.1062, 2.0593, 2.0195, 1.9852, 1.9554, 1.9292, 1.9059, 1.885, 1.8663, 1.8493, 1.8338, 1.8197],
    [44, 4.0617, 3.2093, 2.8165, 2.5837, 2.4271, 2.3133, 2.2262, 2.1572, 2.1009, 2.0539, 2.014, 1.9797, 1.9499, 1.9236, 1.9002, 1.8794, 1.8606, 1.8436, 1.8281, 1.8139],
    [45, 4.0566, 3.2043, 2.8115, 2.5787, 2.4221, 2.3083, 2.2212, 2.1521, 2.0958, 2.0487, 2.0088, 1.9745, 1.9446, 1.9182, 1.8949, 1.874, 1.8551, 1.8381, 1.8226, 1.8084],
    [46, 4.0518, 3.1996, 2.8068, 2.574, 2.4174, 2.3035, 2.2164, 2.1473, 2.0909, 2.0438, 2.0039, 1.9695, 1.9395, 1.9132, 1.8898, 1.8688, 1.85, 1.8329, 1.8173, 1.8031],
    [47, 4.0471, 3.1951, 2.8024, 2.5695, 2.4128, 2.299, 2.2118, 2.1427, 2.0862, 2.0391, 1.9991, 1.9647, 1.9347, 1.9083, 1.8849, 1.8639, 1.845, 1.8279, 1.8123, 1.798],
    [48, 4.0426, 3.1907, 2.7981, 2.5653, 2.4085, 2.2946, 2.2074, 2.1382, 2.0817, 2.0346, 1.9946, 1.9601, 1.9301, 1.9037, 1.8802, 1.8592, 1.8402, 1.8231, 1.8075, 1.7932],
    [49, 4.0384, 3.1866, 2.7939, 2.5611, 2.4044, 2.2904, 2.2032, 2.134, 2.0774, 2.0303, 1.9902, 1.9558, 1.9257, 1.8992, 1.8757, 1.8547, 1.8357, 1.8185, 1.8029, 1.7886],
    [50, 4.0343, 3.1826, 2.79, 2.5572, 2.4004, 2.2864, 2.1992, 2.1299, 2.0734, 2.0261, 1.9861, 1.9515, 1.9214, 1.8949, 1.8714, 1.8503, 1.8313, 1.8141, 1.7985, 1.7841],
    [51, 4.0303, 3.1788, 2.7862, 2.5534, 2.3966, 2.2826, 2.1953, 2.126, 2.0694, 2.0222, 1.982, 1.9475, 1.9174, 1.8908, 1.8673, 1.8462, 1.8272, 1.8099, 1.7942, 1.7798],
    [52, 4.0266, 3.1752, 2.7826, 2.5498, 2.3929, 2.2789, 2.1916, 2.1223, 2.0656, 2.0184, 1.9782, 1.9436, 1.9134, 1.8869, 1.8633, 1.8422, 1.8231, 1.8059, 1.7901, 1.7758],
    [53, 4.023, 3.1716, 2.7791, 2.5463, 2.3894, 2.2754, 2.1881, 2.1187, 2.062, 2.0147, 1.9745, 1.9399, 1.9097, 1.8831, 1.8595, 1.8383, 1.8193, 1.802, 1.7862, 1.7718],
    [54, 4.0196, 3.1683, 2.7757, 2.5429, 2.3861, 2.272, 2.1846, 2.1152, 2.0585, 2.0112, 1.971, 1.9363, 1.9061, 1.8795, 1.8558, 1.8346, 1.8155, 1.7982, 1.7825, 1.768],
    [55, 4.0162, 3.165, 2.7725, 2.5397, 2.3828, 2.2687, 2.1813, 2.1119, 2.0552, 2.0078, 1.9676, 1.9329, 1.9026, 1.876, 1.8523, 1.8311, 1.812, 1.7946, 1.7788, 1.7644],
    [56, 4.0129, 3.1618, 2.7694, 2.5366, 2.3797, 2.2656, 2.1781, 2.1087, 2.0519, 2.0045, 1.9642, 1.9296, 1.8993, 1.8726, 1.8489, 1.8276, 1.8085, 1.7912, 1.7753, 1.7608],
    [57, 4.0099, 3.1589, 2.7665, 2.5336, 2.3767, 2.2625, 2.1751, 2.1056, 2.0488, 2.0014, 1.9611, 1.9264, 1.896, 1.8693, 1.8456, 1.8244, 1.8052, 1.7878, 1.772, 1.7575],
    [58, 4.0069, 3.1559, 2.7635, 2.5307, 2.3738, 2.2596, 2.1721, 2.1026, 2.0458, 1.9983, 1.958, 1.9233, 1.8929, 1.8662, 1.8424, 1.8212, 1.802, 1.7846, 1.7687, 1.7542],
    [59, 4.0039, 3.1531, 2.7608, 2.5279, 2.371, 2.2568, 2.1693, 2.0997, 2.0429, 1.9954, 1.9551, 1.9203, 1.8899, 1.8632, 1.8394, 1.8181, 1.7989, 1.7815, 1.7656, 1.751],
    [60, 4.0012, 3.1504, 2.7581, 2.5252, 2.3683, 2.254, 2.1665, 2.097, 2.0401, 1.9926, 1.9522, 1.9174, 1.887, 1.8603, 1.8364, 1.8151, 1.7959, 1.7784, 1.7625, 1.748],
    [61, 3.9985, 3.1478, 2.7555, 2.5226, 2.3657, 2.2514, 2.1639, 2.0943, 2.0374, 1.9899, 1.9495, 1.9146, 1.8842, 1.8574, 1.8336, 1.8122, 1.793, 1.7755, 1.7596, 1.745],
    [62, 3.9959, 3.1453, 2.753, 2.5201, 2.3631, 2.2489, 2.1613, 2.0917, 2.0348, 1.9872, 1.9468, 1.9119, 1.8815, 1.8547, 1.8308, 1.8095, 1.7902, 1.7727, 1.7568, 1.7422],
    [63, 3.9934, 3.1428, 2.7506, 2.5176, 2.3607, 2.2464, 2.1588, 2.0892, 2.0322, 1.9847, 1.9442, 1.9093, 1.8789, 1.852, 1.8282, 1.8068, 1.7875, 1.77, 1.754, 1.7394],
    [64, 3.9909, 3.1404, 2.7482, 2.5153, 2.3583, 2.244, 2.1564, 2.0868, 2.0298, 1.9822, 1.9417, 1.9068, 1.8763, 1.8495, 1.8256, 1.8042, 1.7849, 1.7673, 1.7514, 1.7368],
    [65, 3.9885, 3.1381, 2.7459, 2.513, 2.356, 2.2417, 2.1541, 2.0844, 2.0274, 1.9798, 1.9393, 1.9044, 1.8739, 1.847, 1.8231, 1.8017, 1.7823, 1.7648, 1.7488, 1.7342],
    [66, 3.9862, 3.1359, 2.7437, 2.5108, 2.3538, 2.2395, 2.1518, 2.0821, 2.0251, 1.9775, 1.937, 1.902, 1.8715, 1.8446, 1.8207, 1.7992, 1.7799, 1.7623, 1.7463, 1.7316],
    [67, 3.9841, 3.1338, 2.7416, 2.5087, 2.3516, 2.2373, 2.1497, 2.0799, 2.0229, 1.9752, 1.9347, 1.8997, 1.8692, 1.8423, 1.8183, 1.7968, 1.7775, 1.7599, 1.7439, 1.7292],
    [68, 3.9819, 3.1317, 2.7395, 2.5066, 2.3496, 2.2352, 2.1475, 2.0778, 2.0207, 1.973, 1.9325, 1.8975, 1.867, 1.84, 1.816, 1.7945, 1.7752, 1.7576, 1.7415, 1.7268],
    [69, 3.9798, 3.1297, 2.7375, 2.5046, 2.3475, 2.2332, 2.1455, 2.0757, 2.0186, 1.9709, 1.9303, 1.8954, 1.8648, 1.8378, 1.8138, 1.7923, 1.7729, 1.7553, 1.7393, 1.7246],
    [70, 3.9778, 3.1277, 2.7355, 2.5027, 2.3456, 2.2312, 2.1435, 2.0737, 2.0166, 1.9689, 1.9283, 1.8932, 1.8627, 1.8357, 1.8117, 1.7902, 1.7707, 1.7531, 1.7371, 1.7223],
    [71, 3.9758, 3.1258, 2.7336, 2.5007, 2.3437, 2.2293, 2.1415, 2.0717, 2.0146, 1.9669, 1.9263, 1.8912, 1.8606, 1.8336, 1.8096, 1.7881, 1.7686, 1.751, 1.7349, 1.7202],
    [72, 3.9739, 3.1239, 2.7318, 2.4989, 2.3418, 2.2274, 2.1397, 2.0698, 2.0127, 1.9649, 1.9243, 1.8892, 1.8586, 1.8316, 1.8076, 1.786, 1.7666, 1.7489, 1.7328, 1.7181],
    [73, 3.9721, 3.1221, 2.73, 2.4971, 2.34, 2.2256, 2.1378, 2.068, 2.0108, 1.9631, 1.9224, 1.8873, 1.8567, 1.8297, 1.8056, 1.784, 1.7646, 1.7469, 1.7308, 1.716],
    [74, 3.9703, 3.1204, 2.7283, 2.4954, 2.3383, 2.2238, 2.1361, 2.0662, 2.009, 1.9612, 1.9205, 1.8854, 1.8548, 1.8278, 1.8037, 1.7821, 1.7626, 1.7449, 1.7288, 1.714],
    [75, 3.9685, 3.1186, 2.7266, 2.4937, 2.3366, 2.2221, 2.1343, 2.0645, 2.0073, 1.9595, 1.9188, 1.8836, 1.853, 1.8259, 1.8018, 1.7802, 1.7607, 1.7431, 1.7269, 1.7121],
    [76, 3.9668, 3.117, 2.7249, 2.4921, 2.3349, 2.2204, 2.1326, 2.0627, 2.0055, 1.9577, 1.917, 1.8819, 1.8512, 1.8241, 1.8, 1.7784, 1.7589, 1.7412, 1.725, 1.7102],
    [77, 3.9651, 3.1154, 2.7233, 2.4904, 2.3333, 2.2188, 2.131, 2.0611, 2.0039, 1.956, 1.9153, 1.8801, 1.8494, 1.8223, 1.7982, 1.7766, 1.7571, 1.7394, 1.7232, 1.7084],
    [78, 3.9635, 3.1138, 2.7218, 2.4889, 2.3318, 2.2172, 2.1294, 2.0595, 2.0022, 1.9544, 1.9136, 1.8785, 1.8478, 1.8206, 1.7965, 1.7749, 1.7554, 1.7376, 1.7214, 1.7066],
    [79, 3.9619, 3.1123, 2.7203, 2.4874, 2.3302, 2.2157, 2.1279, 2.0579, 2.0006, 1.9528, 1.912, 1.8769, 1.8461, 1.819, 1.7948, 1.7732, 1.7537, 1.7359, 1.7197, 1.7048],
    [80, 3.9604, 3.1107, 2.7188, 2.4859, 2.3287, 2.2142, 2.1263, 2.0564, 1.9991, 1.9512, 1.9105, 1.8753, 1.8445, 1.8174, 1.7932, 1.7716, 1.752, 1.7342, 1.718, 1.7032],
    [81, 3.9589, 3.1093, 2.7173, 2.4845, 2.3273, 2.2127, 2.1248, 2.0549, 1.9976, 1.9497, 1.9089, 1.8737, 1.8429, 1.8158, 1.7916, 1.77, 1.7504, 1.7326, 1.7164, 1.7015],
    [82, 3.9574, 3.1079, 2.716, 2.483, 2.3258, 2.2113, 2.1234, 2.0534, 1.9962, 1.9482, 1.9074, 1.8722, 1.8414, 1.8143, 1.7901, 1.7684, 1.7488, 1.731, 1.7148, 1.6999],
    [83, 3.956, 3.1065, 2.7146, 2.4817, 2.3245, 2.2099, 2.122, 2.052, 1.9947, 1.9468, 1.906, 1.8707, 1.8399, 1.8127, 1.7886, 1.7669, 1.7473, 1.7295, 1.7132, 1.6983],
    [84, 3.9546, 3.1051, 2.7132, 2.4803, 2.3231, 2.2086, 2.1206, 2.0506, 1.9933, 1.9454, 1.9045, 1.8693, 1.8385, 1.8113, 1.7871, 1.7654, 1.7458, 1.728, 1.7117, 1.6968],
    [85, 3.9532, 3.1039, 2.7119, 2.479, 2.3218, 2.2072, 2.1193, 2.0493, 1.9919, 1.944, 1.9031, 1.8679, 1.8371, 1.8099, 1.7856, 1.7639, 1.7443, 1.7265, 1.7102, 1.6953],
    [86, 3.9519, 3.1026, 2.7106, 2.4777, 2.3205, 2.2059, 2.118, 2.048, 1.9906, 1.9426, 1.9018, 1.8665, 1.8357, 1.8085, 1.7842, 1.7625, 1.7429, 1.725, 1.7088, 1.6938],
    [87, 3.9506, 3.1013, 2.7094, 2.4765, 2.3193, 2.2047, 2.1167, 2.0467, 1.9893, 1.9413, 1.9005, 1.8652, 1.8343, 1.8071, 1.7829, 1.7611, 1.7415, 1.7236, 1.7073, 1.6924],
    [88, 3.9493, 3.1001, 2.7082, 2.4753, 2.318, 2.2034, 2.1155, 2.0454, 1.9881, 1.94, 1.8992, 1.8639, 1.833, 1.8058, 1.7815, 1.7598, 1.7401, 1.7223, 1.706, 1.691],
    [89, 3.9481, 3.0988, 2.707, 2.4741, 2.3169, 2.2022, 2.1143, 2.0442, 1.9868, 1.9388, 1.8979, 1.8626, 1.8317, 1.8045, 1.7802, 1.7584, 1.7388, 1.7209, 1.7046, 1.6896],
    [90, 3.9469, 3.0977, 2.7058, 2.4729, 2.3157, 2.2011, 2.1131, 2.043, 1.9856, 1.9376, 1.8967, 1.8613, 1.8305, 1.8032, 1.7789, 1.7571, 1.7375, 1.7196, 1.7033, 1.6883],
    [91, 3.9457, 3.0965, 2.7047, 2.4718, 2.3146, 2.1999, 2.1119, 2.0418, 1.9844, 1.9364, 1.8955, 1.8601, 1.8292, 1.802, 1.7777, 1.7559, 1.7362, 1.7183, 1.702, 1.687],
    [92, 3.9446, 3.0955, 2.7036, 2.4707, 2.3134, 2.1988, 2.1108, 2.0407, 1.9833, 1.9352, 1.8943, 1.8589, 1.828, 1.8008, 1.7764, 1.7546, 1.735, 1.717, 1.7007, 1.6857],
    [93, 3.9435, 3.0944, 2.7025, 2.4696, 2.3123, 2.1977, 2.1097, 2.0395, 1.9821, 1.934, 1.8931, 1.8578, 1.8269, 1.7996, 1.7753, 1.7534, 1.7337, 1.7158, 1.6995, 1.6845],
    [94, 3.9423, 3.0933, 2.7014, 2.4685, 2.3113, 2.1966, 2.1086, 2.0385, 1.981, 1.9329, 1.892, 1.8566, 1.8257, 1.7984, 1.7741, 1.7522, 1.7325, 1.7146, 1.6982, 1.6832],
    [95, 3.9412, 3.0922, 2.7004, 2.4675, 2.3102, 2.1955, 2.1075, 2.0374, 1.9799, 1.9318, 1.8909, 1.8555, 1.8246, 1.7973, 1.7729, 1.7511, 1.7314, 1.7134, 1.6971, 1.682],
    [96, 3.9402, 3.0912, 2.6994, 2.4665, 2.3092, 2.1945, 2.1065, 2.0363, 1.9789, 1.9308, 1.8898, 1.8544, 1.8235, 1.7961, 1.7718, 1.75, 1.7302, 1.7123, 1.6959, 1.6809],
    [97, 3.9392, 3.0902, 2.6984, 2.4655, 2.3082, 2.1935, 2.1054, 2.0353, 1.9778, 1.9297, 1.8888, 1.8533, 1.8224, 1.7951, 1.7707, 1.7488, 1.7291, 1.7112, 1.6948, 1.6797],
    [98, 3.9381, 3.0892, 2.6974, 2.4645, 2.3072, 2.1925, 2.1044, 2.0343, 1.9768, 1.9287, 1.8877, 1.8523, 1.8213, 1.794, 1.7696, 1.7478, 1.728, 1.71, 1.6936, 1.6786],
    [99, 3.9371, 3.0882, 2.6965, 2.4636, 2.3062, 2.1916, 2.1035, 2.0333, 1.9758, 1.9277, 1.8867, 1.8513, 1.8203, 1.7929, 1.7686, 1.7467, 1.7269, 1.709, 1.6926, 1.6775],
    [100, 3.9361, 3.0873, 2.6955, 2.4626, 2.3053, 2.1906, 2.1025, 2.0323, 1.9748, 1.9267, 1.8857, 1.8502, 1.8193, 1.7919, 1.7675, 1.7456, 1.7259, 1.7079, 1.6915, 1.6764],
    [101, 3.9352, 3.0864, 2.6946, 2.4617, 2.3044, 2.1897, 2.1016, 2.0314, 1.9739, 1.9257, 1.8847, 1.8493, 1.8183, 1.7909, 1.7665, 1.7446, 1.7248, 1.7069, 1.6904, 1.6754],
    [102, 3.9342, 3.0854, 2.6937, 2.4608, 2.3035, 2.1888, 2.1007, 2.0304, 1.9729, 1.9248, 1.8838, 1.8483, 1.8173, 1.7899, 1.7655, 1.7436, 1.7238, 1.7058, 1.6894, 1.6744],
    [103, 3.9333, 3.0846, 2.6928, 2.4599, 2.3026, 2.1879, 2.0997, 2.0295, 1.972, 1.9238, 1.8828, 1.8474, 1.8163, 1.789, 1.7645, 1.7427, 1.7229, 1.7048, 1.6884, 1.6733],
    [104, 3.9325, 3.0837, 2.692, 2.4591, 2.3017, 2.187, 2.0989, 2.0287, 1.9711, 1.9229, 1.8819, 1.8464, 1.8154, 1.788, 1.7636, 1.7417, 1.7219, 1.7039, 1.6874, 1.6723],
    [105, 3.9316, 3.0828, 2.6912, 2.4582, 2.3009, 2.1861, 2.098, 2.0278, 1.9702, 1.922, 1.881, 1.8455, 1.8145, 1.7871, 1.7627, 1.7407, 1.7209, 1.7029, 1.6865, 1.6714],
    [106, 3.9307, 3.082, 2.6903, 2.4574, 2.3, 2.1853, 2.0971, 2.0269, 1.9694, 1.9212, 1.8801, 1.8446, 1.8136, 1.7862, 1.7618, 1.7398, 1.72, 1.702, 1.6855, 1.6704],
    [107, 3.9299, 3.0812, 2.6895, 2.4566, 2.2992, 2.1845, 2.0963, 2.0261, 1.9685, 1.9203, 1.8792, 1.8438, 1.8127, 1.7853, 1.7608, 1.7389, 1.7191, 1.7011, 1.6846, 1.6695],
    [108, 3.929, 3.0804, 2.6887, 2.4558, 2.2984, 2.1837, 2.0955, 2.0252, 1.9677, 1.9195, 1.8784, 1.8429, 1.8118, 1.7844, 1.7599, 1.738, 1.7182, 1.7001, 1.6837, 1.6685],
    [109, 3.9282, 3.0796, 2.6879, 2.455, 2.2976, 2.1828, 2.0947, 2.0244, 1.9669, 1.9186, 1.8776, 1.8421, 1.811, 1.7835, 1.7591, 1.7371, 1.7173, 1.6992, 1.6828, 1.6676],
    [110, 3.9274, 3.0788, 2.6872, 2.4542, 2.2968, 2.1821, 2.0939, 2.0236, 1.9661, 1.9178, 1.8767, 1.8412, 1.8102, 1.7827, 1.7582, 1.7363, 1.7164, 1.6984, 1.6819, 1.6667],
    [111, 3.9266, 3.0781, 2.6864, 2.4535, 2.2961, 2.1813, 2.0931, 2.0229, 1.9653, 1.917, 1.8759, 1.8404, 1.8093, 1.7819, 1.7574, 1.7354, 1.7156, 1.6975, 1.681, 1.6659],
    [112, 3.9258, 3.0773, 2.6857, 2.4527, 2.2954, 2.1806, 2.0924, 2.0221, 1.9645, 1.9163, 1.8751, 1.8396, 1.8085, 1.7811, 1.7566, 1.7346, 1.7147, 1.6967, 1.6802, 1.665],
    [113, 3.9251, 3.0766, 2.6849, 2.452, 2.2946, 2.1798, 2.0916, 2.0213, 1.9637, 1.9155, 1.8744, 1.8388, 1.8077, 1.7803, 1.7558, 1.7338, 1.7139, 1.6958, 1.6793, 1.6642],
    [114, 3.9243, 3.0758, 2.6842, 2.4513, 2.2939, 2.1791, 2.0909, 2.0206, 1.963, 1.9147, 1.8736, 1.8381, 1.8069, 1.7795, 1.755, 1.733, 1.7131, 1.695, 1.6785, 1.6633],
    [115, 3.9236, 3.0751, 2.6835, 2.4506, 2.2932, 2.1784, 2.0902, 2.0199, 1.9623, 1.914, 1.8729, 1.8373, 1.8062, 1.7787, 1.7542, 1.7322, 1.7123, 1.6942, 1.6777, 1.6625],
    [116, 3.9228, 3.0744, 2.6828, 2.4499, 2.2925, 2.1777, 2.0895, 2.0192, 1.9615, 1.9132, 1.8721, 1.8365, 1.8054, 1.7779, 1.7534, 1.7314, 1.7115, 1.6934, 1.6769, 1.6617],
    [117, 3.9222, 3.0738, 2.6821, 2.4492, 2.2918, 2.177, 2.0888, 2.0185, 1.9608, 1.9125, 1.8714, 1.8358, 1.8047, 1.7772, 1.7527, 1.7307, 1.7108, 1.6927, 1.6761, 1.6609],
    [118, 3.9215, 3.0731, 2.6815, 2.4485, 2.2912, 2.1763, 2.0881, 2.0178, 1.9601, 1.9118, 1.8707, 1.8351, 1.804, 1.7765, 1.752, 1.7299, 1.71, 1.6919, 1.6754, 1.6602],
    [119, 3.9208, 3.0724, 2.6808, 2.4479, 2.2905, 2.1757, 2.0874, 2.0171, 1.9594, 1.9111, 1.87, 1.8344, 1.8032, 1.7757, 1.7512, 1.7292, 1.7093, 1.6912, 1.6746, 1.6594],
    [120, 3.9202, 3.0718, 2.6802, 2.4472, 2.2899, 2.175, 2.0868, 2.0164, 1.9588, 1.9105, 1.8693, 1.8337, 1.8026, 1.775, 1.7505, 1.7285, 1.7085, 1.6904, 1.6739, 1.6587],
    [121, 3.9194, 3.0712, 2.6795, 2.4466, 2.2892, 2.1744, 2.0861, 2.0158, 1.9581, 1.9098, 1.8686, 1.833, 1.8019, 1.7743, 1.7498, 1.7278, 1.7078, 1.6897, 1.6732, 1.6579],
    [122, 3.9188, 3.0705, 2.6789, 2.446, 2.2886, 2.1737, 2.0855, 2.0151, 1.9575, 1.9091, 1.868, 1.8324, 1.8012, 1.7736, 1.7491, 1.727, 1.7071, 1.689, 1.6724, 1.6572],
    [123, 3.9181, 3.0699, 2.6783, 2.4454, 2.288, 2.1731, 2.0849, 2.0145, 1.9568, 1.9085, 1.8673, 1.8317, 1.8005, 1.773, 1.7484, 1.7264, 1.7064, 1.6883, 1.6717, 1.6565],
    [124, 3.9176, 3.0693, 2.6777, 2.4448, 2.2874, 2.1725, 2.0842, 2.0139, 1.9562, 1.9078, 1.8667, 1.831, 1.7999, 1.7723, 1.7478, 1.7257, 1.7058, 1.6876, 1.6711, 1.6558],
    [125, 3.9169, 3.0687, 2.6771, 2.4442, 2.2868, 2.1719, 2.0836, 2.0133, 1.9556, 1.9072, 1.866, 1.8304, 1.7992, 1.7717, 1.7471, 1.725, 1.7051, 1.6869, 1.6704, 1.6551],
    [126, 3.9163, 3.0681, 2.6765, 2.4436, 2.2862, 2.1713, 2.083, 2.0126, 1.955, 1.9066, 1.8654, 1.8298, 1.7986, 1.771, 1.7464, 1.7244, 1.7044, 1.6863, 1.6697, 1.6544],
    [127, 3.9157, 3.0675, 2.6759, 2.443, 2.2856, 2.1707, 2.0824, 2.0121, 1.9544, 1.906, 1.8648, 1.8291, 1.7979, 1.7704, 1.7458, 1.7237, 1.7038, 1.6856, 1.669, 1.6538],
    [128, 3.9151, 3.0669, 2.6754, 2.4424, 2.285, 2.1701, 2.0819, 2.0115, 1.9538, 1.9054, 1.8642, 1.8285, 1.7974, 1.7698, 1.7452, 1.7231, 1.7031, 1.685, 1.6684, 1.6531],
    [129, 3.9145, 3.0664, 2.6749, 2.4419, 2.2845, 2.1696, 2.0813, 2.0109, 1.9532, 1.9048, 1.8636, 1.828, 1.7967, 1.7692, 1.7446, 1.7225, 1.7025, 1.6843, 1.6677, 1.6525],
    [130, 3.914, 3.0659, 2.6743, 2.4414, 2.2839, 2.169, 2.0807, 2.0103, 1.9526, 1.9042, 1.863, 1.8273, 1.7962, 1.7685, 1.744, 1.7219, 1.7019, 1.6837, 1.6671, 1.6519],
    [131, 3.9134, 3.0653, 2.6737, 2.4408, 2.2834, 2.1685, 2.0802, 2.0098, 1.9521, 1.9037, 1.8624, 1.8268, 1.7956, 1.768, 1.7434, 1.7213, 1.7013, 1.6831, 1.6665, 1.6513],
    [132, 3.9129, 3.0648, 2.6732, 2.4403, 2.2829, 2.168, 2.0796, 2.0092, 1.9515, 1.9031, 1.8619, 1.8262, 1.795, 1.7674, 1.7428, 1.7207, 1.7007, 1.6825, 1.6659, 1.6506],
    [133, 3.9123, 3.0642, 2.6727, 2.4398, 2.2823, 2.1674, 2.0791, 2.0087, 1.951, 1.9026, 1.8613, 1.8256, 1.7944, 1.7668, 1.7422, 1.7201, 1.7001, 1.6819, 1.6653, 1.65],
    [134, 3.9118, 3.0637, 2.6722, 2.4392, 2.2818, 2.1669, 2.0786, 2.0082, 1.9504, 1.902, 1.8608, 1.8251, 1.7939, 1.7662, 1.7416, 1.7195, 1.6995, 1.6813, 1.6647, 1.6494],
    [135, 3.9112, 3.0632, 2.6717, 2.4387, 2.2813, 2.1664, 2.0781, 2.0076, 1.9499, 1.9015, 1.8602, 1.8245, 1.7933, 1.7657, 1.7411, 1.719, 1.6989, 1.6808, 1.6641, 1.6488],
    [136, 3.9108, 3.0627, 2.6712, 2.4382, 2.2808, 2.1659, 2.0775, 2.0071, 1.9494, 1.901, 1.8597, 1.824, 1.7928, 1.7651, 1.7405, 1.7184, 1.6984, 1.6802, 1.6635, 1.6483],
    [137, 3.9102, 3.0622, 2.6707, 2.4378, 2.2803, 2.1654, 2.077, 2.0066, 1.9488, 1.9004, 1.8592, 1.8235, 1.7922, 1.7646, 1.74, 1.7178, 1.6978, 1.6796, 1.663, 1.6477],
    [138, 3.9098, 3.0617, 2.6702, 2.4373, 2.2798, 2.1649, 2.0766, 2.0061, 1.9483, 1.8999, 1.8586, 1.823, 1.7917, 1.7641, 1.7394, 1.7173, 1.6973, 1.6791, 1.6624, 1.6471],
    [139, 3.9092, 3.0613, 2.6697, 2.4368, 2.2794, 2.1644, 2.0761, 2.0056, 1.9478, 1.8994, 1.8581, 1.8224, 1.7912, 1.7635, 1.7389, 1.7168, 1.6967, 1.6785, 1.6619, 1.6466],
    [140, 3.9087, 3.0608, 2.6692, 2.4363, 2.2789, 2.1639, 2.0756, 2.0051, 1.9473, 1.8989, 1.8576, 1.8219, 1.7907, 1.763, 1.7384, 1.7162, 1.6962, 1.678, 1.6613, 1.646],
    [141, 3.9083, 3.0603, 2.6688, 2.4359, 2.2784, 2.1634, 2.0751, 2.0046, 1.9469, 1.8984, 1.8571, 1.8214, 1.7901, 1.7625, 1.7379, 1.7157, 1.6957, 1.6775, 1.6608, 1.6455],
    [142, 3.9078, 3.0598, 2.6683, 2.4354, 2.2779, 2.163, 2.0747, 2.0042, 1.9464, 1.8979, 1.8566, 1.8209, 1.7897, 1.762, 1.7374, 1.7152, 1.6952, 1.6769, 1.6603, 1.645],
    [143, 3.9073, 3.0594, 2.6679, 2.435, 2.2775, 2.1625, 2.0742, 2.0037, 1.9459, 1.8975, 1.8562, 1.8204, 1.7892, 1.7615, 1.7368, 1.7147, 1.6946, 1.6764, 1.6598, 1.6444],
    [144, 3.9068, 3.0589, 2.6675, 2.4345, 2.277, 2.1621, 2.0737, 2.0033, 1.9455, 1.897, 1.8557, 1.82, 1.7887, 1.761, 1.7364, 1.7142, 1.6941, 1.6759, 1.6592, 1.6439],
    [145, 3.9064, 3.0585, 2.667, 2.4341, 2.2766, 2.1617, 2.0733, 2.0028, 1.945, 1.8965, 1.8552, 1.8195, 1.7882, 1.7605, 1.7359, 1.7137, 1.6936, 1.6754, 1.6587, 1.6434],
    [146, 3.906, 3.0581, 2.6666, 2.4337, 2.2762, 2.1612, 2.0728, 2.0024, 1.9445, 1.8961, 1.8548, 1.819, 1.7877, 1.7601, 1.7354, 1.7132, 1.6932, 1.6749, 1.6582, 1.6429],
    [147, 3.9055, 3.0576, 2.6662, 2.4332, 2.2758, 2.1608, 2.0724, 2.0019, 1.9441, 1.8956, 1.8543, 1.8186, 1.7873, 1.7596, 1.7349, 1.7127, 1.6927, 1.6744, 1.6578, 1.6424],
    [148, 3.9051, 3.0572, 2.6657, 2.4328, 2.2753, 2.1604, 2.072, 2.0015, 1.9437, 1.8952, 1.8539, 1.8181, 1.7868, 1.7591, 1.7344, 1.7123, 1.6922, 1.6739, 1.6573, 1.6419],
    [149, 3.9046, 3.0568, 2.6653, 2.4324, 2.2749, 2.1599, 2.0716, 2.0011, 1.9432, 1.8947, 1.8534, 1.8177, 1.7864, 1.7587, 1.734, 1.7118, 1.6917, 1.6735, 1.6568, 1.6414],
    [150, 3.9042, 3.0564, 2.6649, 2.4319, 2.2745, 2.1595, 2.0711, 2.0006, 1.9428, 1.8943, 1.853, 1.8172, 1.7859, 1.7582, 1.7335, 1.7113, 1.6913, 1.673, 1.6563, 1.641],
    [151, 3.9038, 3.056, 2.6645, 2.4315, 2.2741, 2.1591, 2.0707, 2.0002, 1.9424, 1.8939, 1.8526, 1.8168, 1.7855, 1.7578, 1.7331, 1.7109, 1.6908, 1.6726, 1.6558, 1.6405],
    [152, 3.9033, 3.0555, 2.6641, 2.4312, 2.2737, 2.1587, 2.0703, 1.9998, 1.942, 1.8935, 1.8521, 1.8163, 1.785, 1.7573, 1.7326, 1.7104, 1.6904, 1.6721, 1.6554, 1.64],
    [153, 3.903, 3.0552, 2.6637, 2.4308, 2.2733, 2.1583, 2.0699, 1.9994, 1.9416, 1.8931, 1.8517, 1.8159, 1.7846, 1.7569, 1.7322, 1.71, 1.6899, 1.6717, 1.6549, 1.6396],
    [154, 3.9026, 3.0548, 2.6634, 2.4304, 2.2729, 2.1579, 2.0695, 1.999, 1.9412, 1.8926, 1.8513, 1.8155, 1.7842, 1.7565, 1.7318, 1.7096, 1.6895, 1.6712, 1.6545, 1.6391],
    [155, 3.9021, 3.0544, 2.6629, 2.43, 2.2725, 2.1575, 2.0691, 1.9986, 1.9407, 1.8923, 1.8509, 1.8151, 1.7838, 1.7561, 1.7314, 1.7091, 1.6891, 1.6708, 1.654, 1.6387],
    [156, 3.9018, 3.054, 2.6626, 2.4296, 2.2722, 2.1571, 2.0687, 1.9982, 1.9403, 1.8918, 1.8505, 1.8147, 1.7834, 1.7557, 1.7309, 1.7087, 1.6886, 1.6703, 1.6536, 1.6383],
    [157, 3.9014, 3.0537, 2.6622, 2.4293, 2.2717, 2.1568, 2.0684, 1.9978, 1.94, 1.8915, 1.8501, 1.8143, 1.7829, 1.7552, 1.7305, 1.7083, 1.6882, 1.6699, 1.6532, 1.6378],
    [158, 3.901, 3.0533, 2.6618, 2.4289, 2.2714, 2.1564, 2.068, 1.9974, 1.9396, 1.8911, 1.8497, 1.8139, 1.7826, 1.7548, 1.7301, 1.7079, 1.6878, 1.6695, 1.6528, 1.6374],
    [159, 3.9006, 3.0529, 2.6615, 2.4285, 2.271, 2.156, 2.0676, 1.997, 1.9392, 1.8907, 1.8493, 1.8135, 1.7822, 1.7544, 1.7297, 1.7075, 1.6874, 1.6691, 1.6524, 1.637],
    [160, 3.9002, 3.0525, 2.6611, 2.4282, 2.2706, 2.1556, 2.0672, 1.9967, 1.9388, 1.8903, 1.8489, 1.8131, 1.7818, 1.754, 1.7293, 1.7071, 1.687, 1.6687, 1.6519, 1.6366],
    [161, 3.8998, 3.0522, 2.6607, 2.4278, 2.2703, 2.1553, 2.0669, 1.9963, 1.9385, 1.8899, 1.8485, 1.8127, 1.7814, 1.7537, 1.7289, 1.7067, 1.6866, 1.6683, 1.6515, 1.6361],
    [162, 3.8995, 3.0518, 2.6604, 2.4275, 2.27, 2.155, 2.0665, 1.9959, 1.9381, 1.8895, 1.8482, 1.8124, 1.781, 1.7533, 1.7285, 1.7063, 1.6862, 1.6679, 1.6511, 1.6357],
    [163, 3.8991, 3.0515, 2.6601, 2.4271, 2.2696, 2.1546, 2.0662, 1.9956, 1.9377, 1.8892, 1.8478, 1.812, 1.7806, 1.7529, 1.7282, 1.7059, 1.6858, 1.6675, 1.6507, 1.6353],
    [164, 3.8987, 3.0512, 2.6597, 2.4268, 2.2693, 2.1542, 2.0658, 1.9953, 1.9374, 1.8888, 1.8474, 1.8116, 1.7803, 1.7525, 1.7278, 1.7055, 1.6854, 1.6671, 1.6503, 1.6349],
    [165, 3.8985, 3.0508, 2.6594, 2.4264, 2.2689, 2.1539, 2.0655, 1.9949, 1.937, 1.8885, 1.8471, 1.8112, 1.7799, 1.7522, 1.7274, 1.7052, 1.685, 1.6667, 1.6499, 1.6345],
    [166, 3.8981, 3.0505, 2.6591, 2.4261, 2.2686, 2.1536, 2.0651, 1.9945, 1.9367, 1.8881, 1.8467, 1.8109, 1.7795, 1.7518, 1.727, 1.7048, 1.6846, 1.6663, 1.6496, 1.6341],
    [167, 3.8977, 3.0502, 2.6587, 2.4258, 2.2683, 2.1533, 2.0648, 1.9942, 1.9363, 1.8878, 1.8464, 1.8105, 1.7792, 1.7514, 1.7266, 1.7044, 1.6843, 1.6659, 1.6492, 1.6338],
    [168, 3.8974, 3.0498, 2.6584, 2.4254, 2.268, 2.1529, 2.0645, 1.9939, 1.936, 1.8874, 1.846, 1.8102, 1.7788, 1.7511, 1.7263, 1.704, 1.6839, 1.6656, 1.6488, 1.6334],
    [169, 3.8971, 3.0495, 2.6581, 2.4251, 2.2676, 2.1526, 2.0641, 1.9936, 1.9357, 1.8871, 1.8457, 1.8099, 1.7785, 1.7507, 1.7259, 1.7037, 1.6835, 1.6652, 1.6484, 1.633],
    [170, 3.8967, 3.0492, 2.6578, 2.4248, 2.2673, 2.1523, 2.0638, 1.9932, 1.9353, 1.8868, 1.8454, 1.8095, 1.7781, 1.7504, 1.7256, 1.7033, 1.6832, 1.6648, 1.6481, 1.6326],
    [171, 3.8965, 3.0488, 2.6575, 2.4245, 2.267, 2.152, 2.0635, 1.9929, 1.935, 1.8864, 1.845, 1.8092, 1.7778, 1.75, 1.7252, 1.703, 1.6828, 1.6645, 1.6477, 1.6323],
    [172, 3.8961, 3.0485, 2.6571, 2.4242, 2.2667, 2.1516, 2.0632, 1.9926, 1.9347, 1.8861, 1.8447, 1.8088, 1.7774, 1.7497, 1.7249, 1.7026, 1.6825, 1.6641, 1.6473, 1.6319],
    [173, 3.8958, 3.0482, 2.6568, 2.4239, 2.2664, 2.1513, 2.0628, 1.9923, 1.9343, 1.8858, 1.8443, 1.8085, 1.7771, 1.7493, 1.7246, 1.7023, 1.6821, 1.6638, 1.647, 1.6316],
    [174, 3.8954, 3.0479, 2.6566, 2.4236, 2.266, 2.151, 2.0626, 1.9919, 1.934, 1.8855, 1.844, 1.8082, 1.7768, 1.749, 1.7242, 1.7019, 1.6818, 1.6634, 1.6466, 1.6312],
    [175, 3.8952, 3.0476, 2.6563, 2.4233, 2.2658, 2.1507, 2.0622, 1.9916, 1.9337, 1.8852, 1.8437, 1.8078, 1.7764, 1.7487, 1.7239, 1.7016, 1.6814, 1.6631, 1.6463, 1.6309],
    [176, 3.8948, 3.0473, 2.6559, 2.423, 2.2655, 2.1504, 2.0619, 1.9913, 1.9334, 1.8848, 1.8434, 1.8075, 1.7761, 1.7483, 1.7236, 1.7013, 1.6811, 1.6628, 1.646, 1.6305],
    [177, 3.8945, 3.047, 2.6556, 2.4227, 2.2652, 2.1501, 2.0616, 1.991, 1.9331, 1.8845, 1.8431, 1.8072, 1.7758, 1.748, 1.7232, 1.7009, 1.6808, 1.6624, 1.6456, 1.6302],
    [178, 3.8943, 3.0467, 2.6554, 2.4224, 2.2649, 2.1498, 2.0613, 1.9907, 1.9328, 1.8842, 1.8428, 1.8069, 1.7755, 1.7477, 1.7229, 1.7006, 1.6805, 1.6621, 1.6453, 1.6298],
    [179, 3.8939, 3.0465, 2.6551, 2.4221, 2.2646, 2.1495, 2.0611, 1.9904, 1.9325, 1.8839, 1.8425, 1.8066, 1.7752, 1.7474, 1.7226, 1.7003, 1.6801, 1.6618, 1.645, 1.6295],
    [180, 3.8936, 3.0462, 2.6548, 2.4218, 2.2643, 2.1492, 2.0608, 1.9901, 1.9322, 1.8836, 1.8422, 1.8063, 1.7749, 1.7471, 1.7223, 1.7, 1.6798, 1.6614, 1.6446, 1.6292],
    [181, 3.8933, 3.0458, 2.6545, 2.4216, 2.264, 2.149, 2.0605, 1.9899, 1.9319, 1.8833, 1.8419, 1.806, 1.7746, 1.7468, 1.7219, 1.6997, 1.6795, 1.6611, 1.6443, 1.6289],
    [182, 3.8931, 3.0456, 2.6543, 2.4213, 2.2638, 2.1487, 2.0602, 1.9896, 1.9316, 1.883, 1.8416, 1.8057, 1.7743, 1.7465, 1.7217, 1.6994, 1.6792, 1.6608, 1.644, 1.6286],
    [183, 3.8928, 3.0453, 2.654, 2.421, 2.2635, 2.1484, 2.0599, 1.9893, 1.9313, 1.8827, 1.8413, 1.8054, 1.774, 1.7462, 1.7214, 1.6991, 1.6789, 1.6605, 1.6437, 1.6282],
    [184, 3.8925, 3.045, 2.6537, 2.4207, 2.2632, 2.1481, 2.0596, 1.989, 1.9311, 1.8825, 1.841, 1.8051, 1.7737, 1.7459, 1.721, 1.6987, 1.6786, 1.6602, 1.6434, 1.6279],
    [185, 3.8923, 3.0448, 2.6534, 2.4205, 2.263, 2.1479, 2.0594, 1.9887, 1.9308, 1.8822, 1.8407, 1.8048, 1.7734, 1.7456, 1.7208, 1.6984, 1.6783, 1.6599, 1.643, 1.6276],
    [186, 3.892, 3.0445, 2.6531, 2.4202, 2.2627, 2.1476, 2.0591, 1.9885, 1.9305, 1.8819, 1.8404, 1.8045, 1.7731, 1.7453, 1.7205, 1.6981, 1.678, 1.6596, 1.6428, 1.6273],
    [187, 3.8917, 3.0442, 2.6529, 2.4199, 2.2624, 2.1473, 2.0588, 1.9882, 1.9302, 1.8816, 1.8401, 1.8042, 1.7728, 1.745, 1.7202, 1.6979, 1.6777, 1.6593, 1.6424, 1.627],
    [188, 3.8914, 3.044, 2.6526, 2.4197, 2.2621, 2.1471, 2.0586, 1.9879, 1.9299, 1.8814, 1.8399, 1.804, 1.7725, 1.7447, 1.7199, 1.6976, 1.6774, 1.659, 1.6421, 1.6267],
    [189, 3.8912, 3.0437, 2.6524, 2.4195, 2.2619, 2.1468, 2.0583, 1.9877, 1.9297, 1.8811, 1.8396, 1.8037, 1.7722, 1.7444, 1.7196, 1.6973, 1.6771, 1.6587, 1.6418, 1.6264],
    [190, 3.8909, 3.0435, 2.6521, 2.4192, 2.2617, 2.1466, 2.0581, 1.9874, 1.9294, 1.8808, 1.8393, 1.8034, 1.772, 1.7441, 1.7193, 1.697, 1.6768, 1.6584, 1.6416, 1.6261],
    [191, 3.8906, 3.0432, 2.6519, 2.4189, 2.2614, 2.1463, 2.0578, 1.9871, 1.9292, 1.8805, 1.8391, 1.8032, 1.7717, 1.7439, 1.719, 1.6967, 1.6765, 1.6581, 1.6413, 1.6258],
    [192, 3.8903, 3.043, 2.6516, 2.4187, 2.2611, 2.1461, 2.0575, 1.9869, 1.9289, 1.8803, 1.8388, 1.8029, 1.7714, 1.7436, 1.7188, 1.6964, 1.6762, 1.6578, 1.641, 1.6255],
    [193, 3.8901, 3.0427, 2.6514, 2.4184, 2.2609, 2.1458, 2.0573, 1.9866, 1.9286, 1.88, 1.8385, 1.8026, 1.7712, 1.7433, 1.7185, 1.6961, 1.6759, 1.6575, 1.6407, 1.6252],
    [194, 3.8899, 3.0425, 2.6512, 2.4182, 2.2606, 2.1456, 2.057, 1.9864, 1.9284, 1.8798, 1.8383, 1.8023, 1.7709, 1.7431, 1.7182, 1.6959, 1.6757, 1.6572, 1.6404, 1.6249],
    [195, 3.8896, 3.0422, 2.6509, 2.418, 2.2604, 2.1453, 2.0568, 1.9861, 1.9281, 1.8795, 1.838, 1.8021, 1.7706, 1.7428, 1.7179, 1.6956, 1.6754, 1.657, 1.6401, 1.6247],
    [196, 3.8893, 3.042, 2.6507, 2.4177, 2.2602, 2.1451, 2.0566, 1.9859, 1.9279, 1.8793, 1.8377, 1.8018, 1.7704, 1.7425, 1.7177, 1.6953, 1.6751, 1.6567, 1.6399, 1.6244],
    [197, 3.8891, 3.0418, 2.6504, 2.4175, 2.26, 2.1448, 2.0563, 1.9856, 1.9277, 1.879, 1.8375, 1.8016, 1.7701, 1.7423, 1.7174, 1.6951, 1.6748, 1.6564, 1.6396, 1.6241],
    [198, 3.8889, 3.0415, 2.6502, 2.4173, 2.2597, 2.1446, 2.0561, 1.9854, 1.9274, 1.8788, 1.8373, 1.8013, 1.7699, 1.742, 1.7172, 1.6948, 1.6746, 1.6562, 1.6393, 1.6238],
    [199, 3.8886, 3.0413, 2.65, 2.417, 2.2595, 2.1444, 2.0558, 1.9852, 1.9272, 1.8785, 1.837, 1.8011, 1.7696, 1.7418, 1.7169, 1.6946, 1.6743, 1.6559, 1.6391, 1.6236],
    [200, 3.8883, 3.041, 2.6497, 2.4168, 2.2592, 2.1441, 2.0556, 1.9849, 1.9269, 1.8783, 1.8368, 1.8008, 1.7694, 1.7415, 1.7166, 1.6943, 1.6741, 1.6557, 1.6388, 1.62]])
        return ftest[row][col]"""


    def magic_read(self,infile):
        """ 
        reads  a Magic template file, puts data in a list of dictionaries
        """
        print "calling magic_read(self, infile)", infile
        hold,magic_data,magic_record,magic_keys=[],[],{},[]
        try:
            f=open(infile,"rU")
        except:
            return [],'bad_file'
        d = f.readline()[:-1].strip('\n')
        if d[0]=="s" or d[1]=="s":
            delim='space'
        elif d[0]=="t" or d[1]=="t":
            delim='tab'
        else: 
            print 'error reading ', infile
            sys.exit()
        if delim=='space':file_type=d.split()[1]
        if delim=='tab':file_type=d.split('\t')[1]
        if file_type=='delimited':
            if delim=='space':file_type=d.split()[2]
            if delim=='tab':file_type=d.split('\t')[2]
        if delim=='space':line =f.readline()[:-1].split()
        if delim=='tab':line =f.readline()[:-1].split('\t')
        for key in line:
            magic_keys.append(key)
        lines=f.readlines()
        for line in lines[:-1]:
            line.replace('\n','')
            if delim=='space':rec=line[:-1].split()
            if delim=='tab':rec=line[:-1].split('\t')
            hold.append(rec)
        line = lines[-1].replace('\n','')
        if delim=='space':rec=line[:-1].split()
        if delim=='tab':rec=line.split('\t')
        hold.append(rec)
        for rec in hold:
            magic_record={}
            if len(magic_keys) != len(rec):
                
                print "Warning: Uneven record lengths detected: "
                print magic_keys
                print rec
            for k in range(len(rec)):
               magic_record[magic_keys[k]]=rec[k].strip('\n')
            magic_data.append(magic_record)
        magictype=file_type.lower().split("_")
        Types=['er','magic','pmag','rmag']
        if magictype in Types:file_type=file_type.lower()
        print "magic data:"
        print str(magic_data)[:500] + "..."
        print "file_type", file_type
        return magic_data,file_type
    

    def get_specs(self,data):
        """
         takes a magic format file and returns a list of unique specimen names
        """
    # sort the specimen names
    #
        print "calling get_specs()"
        speclist=[]
        for rec in data:
          spec=rec["er_specimen_name"]
          if spec not in speclist:speclist.append(spec)
        speclist.sort()
        return speclist
    


    def sortarai(self,datablock,s,Zdiff):
        """
         sorts data block in to first_Z, first_I, etc.
        """
        print "calling sortarai()"
        first_Z,first_I,zptrm_check,ptrm_check,ptrm_tail=[],[],[],[],[]
        field,phi,theta="","",""
        starthere=0
        Treat_I,Treat_Z,Treat_PZ,Treat_PI,Treat_M=[],[],[],[],[]
        ISteps,ZSteps,PISteps,PZSteps,MSteps=[],[],[],[],[]
        GammaChecks=[] # comparison of pTRM direction acquired and lab field
        Mkeys=['measurement_magn_moment','measurement_magn_volume','measurement_magn_mass','measurement_magnitude']
        rec=datablock[0]
        for key in Mkeys:
            if key in rec.keys() and rec[key]!="":
                momkey=key
                break
    # first find all the steps
        for k in range(len(datablock)):
            rec=datablock[k]
            if "treatment_temp" in rec.keys():
                temp=float(rec["treatment_temp"])
            elif "treatment_mw_power" in rec.keys():
                temp=float(rec["treatment_mw_power"])
                
            methcodes=[]
            tmp=rec["magic_method_codes"].split(":")
            for meth in tmp:
                methcodes.append(meth.strip())
            # for thellier-thellier
            if 'LT-T-I' in methcodes and 'LP-PI-TRM' in methcodes and 'LP-TRM' not in methcodes :
                Treat_I.append(temp)
                ISteps.append(k)
                if field=="":field=float(rec["treatment_dc_field"])
                if phi=="":
                    phi=float(rec['treatment_dc_field_phi'])
                    theta=float(rec['treatment_dc_field_theta'])
                    
            # for Microwave
            if 'LT-M-I' in methcodes and 'LP-PI-M' in methcodes :
                Treat_I.append(temp)
                ISteps.append(k)
                if field=="":field=float(rec["treatment_dc_field"])
                if phi=="":
                    phi=float(rec['treatment_dc_field_phi'])
                    theta=float(rec['treatment_dc_field_theta'])

    # stick  first zero field stuff into first_Z 
            if 'LT-NO' in methcodes:
                Treat_Z.append(temp)
                ZSteps.append(k)
            if 'LT-T-Z' in methcodes or 'LT-M-Z' in methcodes: 
                Treat_Z.append(temp)
                ZSteps.append(k)
            if 'LT-PTRM-Z' :
                Treat_PZ.append(temp)
                PZSteps.append(k)
            if 'LT-PTRM-I' in methcodes or 'LT-PMRM-I' in methcodes:
                Treat_PI.append(temp)
                PISteps.append(k)
            if 'LT-PTRM-MD' in methcodes:
                Treat_M.append(temp)
                MSteps.append(k)
            if 'LT-NO' in methcodes:
                dec=float(rec["measurement_dec"])
                inc=float(rec["measurement_inc"])
                str=float(rec[momkey])
                if 'LP-PI-M'  not in methcodes:
                    first_I.append([273,0.,0.,0.,1])
                    first_Z.append([273,dec,inc,str,1])  # NRM step
                else:
                    first_I.append([0,0.,0.,0.,1])
                    first_Z.append([0,dec,inc,str,1])  # NRM step
                    
        #---------------------
        # find  IZ and ZI
        #---------------------
                    
                
        for temp in Treat_I: # look through infield steps and find matching Z step
            if temp in Treat_Z: # found a match
                istep=ISteps[Treat_I.index(temp)]
                irec=datablock[istep]
                methcodes=[]
                tmp=irec["magic_method_codes"].split(":")
                for meth in tmp: methcodes.append(meth.strip())
                brec=datablock[istep-1] # take last record as baseline to subtract  
                zstep=ZSteps[Treat_Z.index(temp)]
                zrec=datablock[zstep]
        # sort out first_Z records 
                if "LP-PI-TRM-IZ" in methcodes or "LP-PI-M-IZ" in methcodes: 
                    ZI=0    
                else:   
                    ZI=1    
                dec=float(zrec["measurement_dec"])
                inc=float(zrec["measurement_inc"])
                str=float(zrec[momkey])
                first_Z.append([temp,dec,inc,str,ZI])
        # sort out first_I records 
                idec=float(irec["measurement_dec"])
                iinc=float(irec["measurement_inc"])
                istr=float(irec[momkey])
                X=self.dir2cart([idec,iinc,istr])
                BL=self.dir2cart([dec,inc,str])
                I=[]
                for c in range(3): I.append((X[c]-BL[c]))
                if I[2]!=0:
                    iDir=self.cart2dir(I)
                    if Zdiff==0:
                        first_I.append([temp,iDir[0],iDir[1],iDir[2],ZI])
                    else:
                        first_I.append([temp,0.,0.,I[2],ZI])
##                    gamma=angle([iDir[0],iDir[1]],[phi,theta])
                else:
                    first_I.append([temp,0.,0.,0.,ZI])
##                    gamma=0.0
##    # put in Gamma check (infield trm versus lab field)
##                if 180.-gamma<gamma:
##                    gamma=180.-gamma
##                GammaChecks.append([temp-273.,gamma])


        #---------------------
        # find Thellier Thellier protocol
        #---------------------
        if 'LP-PI-II'in methcodes or 'LP-PI-T-II' in methcodes or 'LP-PI-M-II' in methcodes:
            for i in range(1,len(Treat_I)): # look through infield steps and find matching Z step
                if Treat_I[i] == Treat_I[i-1]:
                    # ignore, if there are more than 
                    temp= Treat_I[i]
                    irec1=datablock[ISteps[i-1]]
                    dec1=float(irec1["measurement_dec"])
                    inc1=float(irec1["measurement_inc"])
                    moment1=float(irec1["measurement_magn_moment"])
                    if len(first_I)<2:
                        dec_initial=dec1;inc_initial=inc1
                    cart1=array(self.dir2cart([dec1,inc1,moment1]))
                    irec2=datablock[ISteps[i]]
                    dec2=float(irec2["measurement_dec"])
                    inc2=float(irec2["measurement_inc"])
                    moment2=float(irec2["measurement_magn_moment"])
                    cart2=array(self.dir2cart([dec2,inc2,moment2]))

                    # check if its in the same treatment
                    if Treat_I[i] == Treat_I[i-2] and dec2!=dec_initial and inc2!=inc_initial:
                        continue
                    if dec1!=dec2 and inc1!=inc2:
                        zerofield=(cart2+cart1)/2
                        infield=(cart2-cart1)/2

                        DIR_zerofield=self.cart2dir(zerofield)
                        DIR_infield=self.cart2dir(infield)

                        first_Z.append([temp,DIR_zerofield[0],DIR_zerofield[1],DIR_zerofield[2],0])
                        first_I.append([temp,DIR_infield[0],DIR_infield[1],DIR_infield[2],0])
    

        #---------------------
        # find  pTRM checks
        #---------------------
                    
        for temp in Treat_PI: # look through infield steps and find matching Z step
            if 'LP-PI-II' not in methcodes:
                step=PISteps[Treat_PI.index(temp)]
                rec=datablock[step]
                dec=float(rec["measurement_dec"])
                inc=float(rec["measurement_inc"])
                str=float(rec[momkey])
                brec=datablock[step-1] # take last record as baseline to subtract
                pdec=float(brec["measurement_dec"])
                pinc=float(brec["measurement_inc"])
                pint=float(brec[momkey])
                X=self.dir2cart([dec,inc,str])
                prevX=self.dir2cart([pdec,pinc,pint])
                I=[]
                for c in range(3): I.append(X[c]-prevX[c])
                dir1=self.cart2dir(I)
                if Zdiff==0:
                    ptrm_check.append([temp,dir1[0],dir1[1],dir1[2]])
                else:
                    ptrm_check.append([temp,0.,0.,I[2]])
            else:
                step=PISteps[Treat_PI.index(temp)]
                rec=datablock[step]
                dec=float(rec["measurement_dec"])
                inc=float(rec["measurement_inc"])
                moment=float(rec["measurement_magn_moment"])
                for zerofield in first_Z:
                    if zerofield[0]==temp:
                        M1=array(self.dir2cart([dec,inc,moment]))
                        M2=array(self.dir2cart([zerofield[1],zerofield[2],zerofield[3]]))
                        diff=M1-M2
                        diff_cart=self.cart2dir(diff)
                        ptrm_check.append([temp,diff_cart[0],diff_cart[1],diff_cart[2]])
                        
                        
                        
    # in case there are zero-field pTRM checks (not the SIO way)
        for temp in Treat_PZ:
            step=PZSteps[Treat_PZ.index(temp)]
            rec=datablock[step]
            dec=float(rec["measurement_dec"])
            inc=float(rec["measurement_inc"])
            str=float(rec[momkey])
            brec=datablock[step-1]
            pdec=float(brec["measurement_dec"])
            pinc=float(brec["measurement_inc"])
            pint=float(brec[momkey])
            X=self.dir2cart([dec,inc,str])
            prevX=self.dir2cart([pdec,pinc,pint])
            I=[]
            for c in range(3): I.append(X[c]-prevX[c])
            dir2=self.cart2dir(I)
            zptrm_check.append([temp,dir2[0],dir2[1],dir2[2]])
        ## get pTRM tail checks together -
        for temp in Treat_M:
            step=MSteps[Treat_M.index(temp)] # tail check step - just do a difference in magnitude!
            rec=datablock[step]
            str=float(rec[momkey])
            if temp in Treat_Z:
                step=ZSteps[Treat_Z.index(temp)]
                brec=datablock[step]
                pint=float(brec[momkey])
                ptrm_tail.append([temp,0,0,str-pint])  # difference - if negative, negative tail!
            else:
                print s, '  has a tail check with no first zero field step - check input file! for step',temp-273.
    #
    # final check
    #
        if len(first_Z)!=len(first_I):
                   print len(first_Z),len(first_I)
                   print " Something wrong with this specimen! Better fix it or delete it "
                   raw_input(" press return to acknowledge message")
        araiblock=(first_Z,first_I,ptrm_check,ptrm_tail,zptrm_check,GammaChecks)
    
        return araiblock,field


#if __name__ == '__main__':
gui = Arai_GUI()
specimens = gui.Data.keys()
import spd
print specimens
#thing1_tmax = gui.Data[specimens[0]]['t_Arai'][0]
#thing1_tmin = gui.Data[specimens[-1]]['t_Arai'][-1]
#print thing1_tmax, thing1_tmin
#thing1 = spd.PintPars(gui.Data, gui.specimens[0], thing1_tmin, thing1_tmax)
things = []
for n, s in enumerate(specimens):
    print "looping: "
    print s
    print gui.Data[s]['t_Arai']
    tmin = gui.Data[s]['t_Arai'][0]
    tmax = gui.Data[s]['t_Arai'][-1]
    print "tmin is: %s" %(tmin)
    print "tmax is: %s" %(tmax)
    thing = spd.PintPars(gui.Data, s, tmin, tmax)
    things.append(thing)

print things
thing = things[0]
thing1 = things[1]
thing2 = things[2]
thing3 = things[3]
thing4 = things[4]
thing5 = things[5]
