    
    
    
    

double interp_eloss(int currindex, double currE){
    differnce = energy[currindex]-energy[currindex-1];
    low_diff = currE - energy[currindex-1];
    high_diff = energy[currindex] - currE;
    high_weight = low_diff/differnce;
    low_weight = high_diff/differnce;
    interp_energy = low_weight*energy[currindex-1]+high_weight*energy[currindex];
    return interp_energy;
}
