#!/usr/bin/env python

#compare.txt: -*- Text -*-  DESCRIPTIVE TEXT.

#Copyright (c) 2014 Lori Jonestrask
#Author: Lori Jonestrask (mintblue87@gmail.com) .



#seg_min=start_pt;
#seg_max=end_pt;
#seg=(seg_min:seg_max);


#Params.n=length(seg);
#Params.Seg_Ends=[seg_min, seg_max];

#% Max and min temperatures
#Params.Tmin=TRMvec(seg_min,1);
#Params.Tmax=TRMvec(seg_max,1);

#X_seg=Params.Xpts(seg);
#Y_seg=Params.Ypts(seg);

#Params.X_seg = X_seg;
#Params.Y_seg = Y_seg;

#Params.xbar=mean(X_seg);
#Params.ybar=mean(Y_seg);
#U=detrend(X_seg,0); % (Xi-Xbar)
#V=detrend(Y_seg,0); % (Yi-Ybar)

#Params.Y_int=mean(Y_seg)-Params.b*mean(X_seg);
#Params.X_int=-Params.Y_int/Params.b;

#Rev_x=(Params.Ypts-Params.Y_int)./Params.b; % The points reflected about the bes-fit line
#Rev_y=Params.b.*Params.Xpts+Params.Y_int;
#Px=(Params.Xpts+Rev_x)./2; % Average the both sets to get the projected points
#Py=(Params.Ypts+Rev_y)./2;


def York_Regression(x_segment, y_segment, x_mean, y_mean, n, lab_dc_field, steps_Arai):
    x_err = x_segment - x_mean
    y_err = y_segment - y_mean
    york_b = -1* sqrt( sum(y_err**2) / sum(x_err**2) )  # averaged slope                                                  
    york_sigma= sqrt ( (2 * sum(y_err**2) - 2*york_b* sum(x_err*y_err)) / ( (n-2) * sum(x_err**2) ) )
    beta_Coe=abs(york_sigma/york_b
    # y_T is the intercept of the extrepolated line                                                                              
    # through the center of mass (see figure 7 in Coe (1978))                                                                     
    y_T = y_mean - (york_b* x_mean)
    x_T = (-1 * y_T) / york_b  # x intercept                                                                                      
    # # calculate the extarplated data points for f and fvds                                                                      
    x_tag=(y_segment - y_T ) / york_b # returns array of y points minus the y intercept, divided by slope                    
    y_tag=york_b*x_segment + y_T
    # intersect of the dashed square and the horizontal dahed line  next to delta-y-5 in figure 7, Coe (1987)
    x_prime=(x_segment+x_tag) / 2
    y_prime=(y_segment+y_tag) / 2
    delta_x_prime = abs(x_prime[-1] - x_prime[0]) #Lj add.  this is the TRM length of the best fit line
    delta_y_prime = abs(y_prime[-1] - y_prime[0]) # LJ add.  this is the NRM length of the best fit line              
    f_Coe = delta_y_prime / y_T  # LJ added



# python   ,    MATLAB
#x_segment,   x_seg
#x_mean    ,   xbar
# x_err     , U  # confirmed, these are the same

# x_tag    , Rev_x  # BUT mine uses y_segment, his uses all ypts
# x_prime    ,  Px # But depends on Rev_x/x_tag, so not same
