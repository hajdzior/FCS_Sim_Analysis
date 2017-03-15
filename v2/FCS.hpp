#ifndef GUARD_system
#define GUARD_system

#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

struct Coordinates
{
    double x;
    double y;
    double z;
};

struct Box
{
    double x_1, x_2, y_1, y_2, z_1, z_2; 

   
    Box (double _x_1, double _x_2, double _y_1, double _y_2, double _z_1, double _z_2) {
        x_1 = _x_1;
        x_2 = _x_2;
        y_1 = _y_1;
        y_2 = _y_2;
        z_1 = _z_1;
        z_2 = _z_2;
    }

    double volume(){
        // Box volume in um^3.
        return (x_2 - x_1) * (y_2 - y_1) * (z_2 - z_1);
    }

    double volume_L(){
        // Box volume in liters.
        return volume() * 1e-15;
    }
};


class GaussianMDF
{
    // This class implements a Gaussian-shaped PSF function
        double w_xy, w_z;

        public:  
        GaussianMDF (double _w_xy, double _w_z) {
            w_xy =  _w_xy;
            w_z =  _w_z;
        }

        double EID(double x, double y, double z) {
            return exp( -(2.0*x*x/(w_xy*w_xy) + 2.0*y*y/(w_xy*w_xy) + 2.0*z*z/(w_z*w_z)) );
        }

        double eval(double x, double y, double z) {
            return exp( -(2.0*x*x/(w_xy*w_xy) + 2.0*y*y/(w_xy*w_xy) + 2.0*z*z/(w_z*w_z)) );
        }
};


class SemigeometricMDF {
    /* This class implements the Molecular Detection Function (MDF)
    within a semigeometric approximation.

    NA:  objective's numerical aperture
    M:  magnification
    n:  refractive index of the immersion medium
    r:  pinhole radius [um]
    lambda_ex:  exciatation wavelength [um]
    lambda_em:  emission
    w_0:  laser beam waist (e-2) radius [um]

    CEF: Collection Effciency Function
    EID: Excitation Intensity Distribution
    */

    double NA;
    double M;
    double n;
    double r;
    double lambda_ex;
    double lambda_em;
    double w_0;
    double phi_f;

    double a;
    double Theta;
    double omega;

public:
    SemigeometricMDF (
        double _NA,
        double _M,
        double _n,
        double _r,
        double _lambda_ex,
        double _lambda_em,
        double _w_0,
        double _phi_f)
    {
        NA = _NA;
        M = _M;
        n = _n;
        r = _r;
        lambda_ex = _lambda_ex;
        lambda_em = _lambda_em;
        w_0 = _w_0;
        phi_f = _phi_f;

        a = r/M;
        Theta = asin(NA/n);  // maximum angle of light collection
        omega = lambda_em/(M_PI * tan(Theta));
    } ;


    SemigeometricMDF ( ) :
        NA(1.2),
        M(60.0),
        n(1.4),
        r(50.0),
        lambda_ex(0.635),
        lambda_em(0.635),
        w_0(0.300),
        phi_f(1.00)
    {
        a = r/M;
        Theta = asin(NA/n);  // maximum angle of light collection
        omega = lambda_em/(M_PI * tan(Theta));
    } ;


    double CEF(double x, double y, double z)
    {
        /* light-Collection Efficiency Function */

        double rho = sqrt(x*x + y*y);
        double R = omega * sqrt( 1.0 + ((z)*lambda_em/(M_PI*omega*omega)) * ((z)*lambda_em/(M_PI*omega*omega)) );
        double rval = 0.0;

        if (rho >= (R + a)) {
            rval = 0.0;
        } else if ( (fabs(R-a) < rho) && (rho < (R+a) ) ) {
            double Theta_1 = acos( (a*a + rho*rho -R*R)/(2.0*a*rho) );
            double Theta_2 = acos( (R*R + rho*rho -a*a)/(2.0*R*rho) );
            double Delta = sqrt( (a+rho+R)*(-a+rho+R)*(a-rho+R)*(a+rho-R) )/2.0;
            rval = std::max(a, omega)*std::max(a, omega) *
                   ( (Theta_1*a*a) + (Theta_2*R*R) - Delta)/(M_PI* (a*a) * (R*R));
        } else if (rho <= fabs(R-a)) {
            rval = std::max(a, omega)*std::max(a, omega)/((std::max(a, R)*std::max(a, R)));
        }
        return rval;
    }


    double EID(double x, double y, double z)
    {
        /* Excitation Intensity Distribution */
        double rho = sqrt( x*x + y*y );
        double w = w_0*sqrt( 1.0 + ( z*lambda_ex/( M_PI*w_0*w_0))*(z*lambda_ex/( M_PI*w_0*w_0 )) );
        // Normalized to 1 at the centre 
        return w_0*w_0 * ( 1.0/(w*w))* exp(-2.0*rho*rho/(w*w) );
    }


    double eval(double x, double y, double z)
    {
        // double k_ex = EID(x, y, z);  // position-dependent maximum excitation rate
        //return phi_f*k_ex*EID(x, y, z)/(1.0 + k_ex*EID(x, y, z)/k_sat)*CEF(x, y, z);
        // dbg AM no saturation effects k_sat =  infinity
        // return phi_f*k_ex*CEF(x, y, z);
        return EID(x, y, z)*CEF(x, y, z);
    }
};

//-----------------------------------------------------------------------------
class MDF
{
public:
    MDF (
        std::string _confocal_geometry,       
        double _w_xy,
        double _w_z,        
        double _NA,
        double _M,
        double _n,
        double _r,
        double _lambda_ex,
        double _lambda_em,
        double _w_0,
        double _phi_f)
    {

        confocal_geometry = _confocal_geometry;
        w_xy = _w_xy;
        w_z = _w_z;       
        NA = _NA;
        M = _M;
        n = _n;
        r = _r;
        lambda_ex = _lambda_ex;
        lambda_em = _lambda_em;
        w_0 = _w_0;
        phi_f = _phi_f;

        a = r/M;
        Theta = asin(NA/n);  // maximum angle of light collection
        omega = lambda_em/(M_PI * tan(Theta));
    } ;

    MDF (
        double _w_xy,
        double _w_z)
    {
        confocal_geometry = "Gaussian";
        w_xy = _w_xy;
        w_z = _w_z;       
    } ;

    MDF (
        double _NA,
        double _M,
        double _n,
        double _r,
        double _lambda_ex,
        double _lambda_em,
        double _w_0,
        double _phi_f)
    {
        confocal_geometry = "Semigeometric";
        NA = _NA;
        M = _M;
        n = _n;
        r = _r;
        lambda_ex = _lambda_ex;
        lambda_em = _lambda_em;
        w_0 = _w_0;
        phi_f = _phi_f;

        a = r/M;
        Theta = asin(NA/n);  // maximum angle of light collection
        omega = lambda_em/(M_PI * tan(Theta));
    } ;

/*    MDF ( ) :*/
/*        type ("Gaussian"),*/
/*        w_xy (0.2),*/
/*        w_z (1.0),        */
/*        NA(1.2),*/
/*        M(60.0),*/
/*        n(1.4),*/
/*        r(50.0),*/
/*        lambda_ex(0.635),*/
/*        lambda_em(0.635),*/
/*        w_0(0.300),*/
/*        phi_f(1.00)*/
/*    {*/
/*        a = r/M;*/
/*        Theta = asin(NA/n);  // maximum angle of light collection*/
/*        omega = lambda_em/(M_PI * tan(Theta));*/
/*    } ;*/

    double EID(double x, double y, double z)
    {   
        double ret = 0.0;
        if (confocal_geometry == "\"Gaussian\"") {
            ret = GaussianEID(x,y,z);           
        }
        else if (confocal_geometry == "\"Semigeometric\"") {  
            ret = SemigeometricEID(x,y,z);    
        }
        else 
        {
            std::cout << "FCS.hpp: double EID: Unknown geometry. Exiting ..." << std::endl; 
        }
        return ret;
    }

    double eval(double x, double y, double z)
    {   
        double ret = 0.0;
        if (confocal_geometry == "\"Gaussian\"") {
            ret = GaussianEval(x,y,z);           
        }
        else if (confocal_geometry == "\"Semigeometric\"") {  
            ret = SemigeometricEval(x,y,z);    
        }
        else 
        {
            std::cout << "FCS.hpp: double eval: Unknown geometry. Exiting ..." << std::endl; 
        }
        return ret;    
    }

private: 
    std::string confocal_geometry;
    double w_xy, w_z;
    double NA;
    double M;
    double n;
    double r;
    double lambda_ex;
    double lambda_em;
    double w_0;
    double phi_f;

    double a;
    double Theta;
    double omega;


    double GaussianEID(double x, double y, double z) {
        return exp( -(2.0*x*x/(w_xy*w_xy) + 2.0*y*y/(w_xy*w_xy) + 2.0*z*z/(w_z*w_z)) );
    }

    double GaussianEval(double x, double y, double z) {
        return exp( -(2.0*x*x/(w_xy*w_xy) + 2.0*y*y/(w_xy*w_xy) + 2.0*z*z/(w_z*w_z)) );
    }


    double SemigeometricCEF(double x, double y, double z)
    {
        /* light Collection Efficiency Function */

        double rho = sqrt(x*x + y*y);
        double R = omega * sqrt( 1.0 + ((z)*lambda_em/(M_PI*omega*omega)) * ((z)*lambda_em/(M_PI*omega*omega)) );
        double rval = 0.0;

        if (rho >= (R + a)) {
            rval = 0.0;
        } else if ( (fabs(R-a) < rho) && (rho < (R+a) ) ) {
            double Theta_1 = acos( (a*a + rho*rho -R*R)/(2.0*a*rho) );
            double Theta_2 = acos( (R*R + rho*rho -a*a)/(2.0*R*rho) );
            double Delta = sqrt( (a+rho+R)*(-a+rho+R)*(a-rho+R)*(a+rho-R) )/2.0;
            rval = std::max(a, omega)*std::max(a, omega) *
                   ( (Theta_1*a*a) + (Theta_2*R*R) - Delta)/(M_PI* (a*a) * (R*R));
        } else if (rho <= fabs(R-a)) {
            rval = std::max(a, omega)*std::max(a, omega)/((std::max(a, R)*std::max(a, R)));
        }
        return rval;
    }

    double SemigeometricEID(double x, double y, double z)
    {
        /* Excitation Intensity Distribution */
        double rho = sqrt( x*x + y*y );
        double w = w_0*sqrt( 1.0 + ( z*lambda_ex/( M_PI*w_0*w_0))*(z*lambda_ex/( M_PI*w_0*w_0 )) );
        // Normalized to 1 at the centre 
        return w_0*w_0 * ( 1.0/(w*w))* exp(-2.0*rho*rho/(w*w) );
    }


    double SemigeometricEval(double x, double y, double z)
    {
        // double k_ex = EID(x, y, z);  // position-dependent maximum excitation rate
        //return phi_f*k_ex*EID(x, y, z)/(1.0 + k_ex*EID(x, y, z)/k_sat)*CEF(x, y, z);
        // dbg AM no saturation effects k_sat =  infinity
        // return phi_f*k_ex*CEF(x, y, z);
        return SemigeometricEID(x, y, z)*SemigeometricCEF(x, y, z);
    }
};


struct Molecules 
{
    double D;
    std::vector<Coordinates> positions;
    std::vector<unsigned int> states;
    std::string blinking;   

    Molecules (double _D, size_t n_particles, Box box, gsl_rng * rng, std::string _blinking) {
        D = _D;
        Coordinates coordinates;
        blinking = _blinking; 
        
        for( size_t i=0; i<n_particles; ++i ) {
            coordinates.x = gsl_rng_uniform(rng)*(box.x_2-box.x_1) + box.x_1;
            coordinates.y = gsl_rng_uniform(rng)*(box.y_2-box.y_1) + box.y_1;
            coordinates.z = gsl_rng_uniform(rng)*(box.z_2-box.z_1) + box.z_1;                       
            positions.push_back(coordinates);
            states.push_back(0);     
        } 
    }
    
    void update_positions_periodic (size_t n_particles, Box box, gsl_rng * rng, double dt) {
        double sigma = sqrt( 2.0 * D * dt );

        for(size_t i=0; i<n_particles; ++i){
                positions[i].x = positions[i].x + gsl_ran_gaussian (rng, sigma);
                positions[i].y = positions[i].y + gsl_ran_gaussian (rng, sigma);
                positions[i].z = positions[i].z + gsl_ran_gaussian (rng, sigma);                          

                // periodic boundary conditions    
                if (positions[i].x < box.x_1) {
                    positions[i].x = positions[i].x + (box.x_2 - box.x_1);
                    states[i] = 0; 
                }
                if (positions[i].x > box.x_2) {
                    positions[i].x = positions[i].x - (box.x_2 - box.x_1);
                    states[i] = 0;
                }
                if (positions[i].y < box.y_1) { 
                    positions[i].y = positions[i].y + (box.y_2 - box.y_1); 
                    states[i] = 0;
                }                    
                if (positions[i].y > box.y_2) {
                    positions[i].y = positions[i].y - (box.y_2 - box.y_1);
                    states[i] = 0;
                }
                if (positions[i].z < box.z_1) { 
                    positions[i].z = positions[i].z + (box.z_2 - box.z_1); 
                    states[i] = 0;
                }                    
                if (positions[i].z > box.z_2) {
                    positions[i].z = positions[i].z - (box.z_2 - box.z_1);
                    states[i] = 0;
                }     
        }
    }

    void update_states (size_t n_particles, 
                        gsl_rng * rng, 
                        double k_ex_max, 
                        double phi_st, 
                        double phi_sb, 
                        double k_ts, 
                        double k_tb, 
                        double dt, 
                        MDF mdf) 
    {
        if (blinking == "\"triplet\"") 
        {        
            // Delare working variables
            double k_st = 0.0;        
            double k_sb = 0.0; 
            double current_ex = 0.0;

            for(size_t i=0; i<n_particles; ++i) {

                assert (k_ex_max*dt < 0.2);
                assert (phi_st + phi_sb <= 1.0);

                current_ex = mdf.EID(positions[i].x, positions[i].y, positions[i].z);
                // position dependent singlet -> triplet rate                     
                k_st = k_ex_max*phi_st*current_ex;
                // singlet -> bleached rate                
                k_sb = k_ex_max*phi_sb*current_ex;
                // triplet -> singlet rate 

                double tmp = gsl_rng_uniform( rng );                        
                if (states[i] == 0) { // singlet
                    if ( tmp > (k_st+k_sb)*dt ) {  
                        continue;                               
                    }                    
                    else if ( tmp > k_sb*dt ) {
                        //std::cout << "s->t " << std::endl;                                                                      
                        states[i] = 1; // triplet   
                    }
                    else {
                        //std::cout << "s->b " << std::endl;                               
                        states[i] = 2; // bleached 
                    }
                }
                else if (states[i] == 1) { // triplet
                    if ( tmp > ( k_ts+k_tb)*dt ) {                               
                        continue;
                    }
                    else if (tmp > k_tb*dt) {
                        //std::cout << "t->s " << std::endl;
                        states[i] = 0; // singlet
                    } 
                    else {
                        //std::cout << "t->b " << std::endl; 
                        states[i] = 2 ;// bleached
                    }
                }
                else if (states[i] == 2) { // bleached
                    continue;
                }
                else {
                    std::cout << "FCS.hpp: state = 3, an error" << std::endl;
                    exit(1);
                }
            }
        }
        else 
        {
            std::cout << "FCS.hpp: Powerlaw blinking is not implemented yet" << std::endl;
            exit(1);
        }
    }

    unsigned int arrivals( size_t n_particles, 
                            double k_ex_max, 
                            double phi_f, 
                            double phi_det, 
                            double dt, 
                            //GaussianMDF mdf,
                            MDF mdf, 
                            gsl_rng * rng ) {
        double k_det = 0.0;  
        unsigned int arrivals = 0;
          
        for(size_t i=0; i<n_particles; ++i) {
            if (states[i] == 0) {
                k_det += mdf.eval(positions[i].x, positions[i].y, positions[i].z);
            }
            //std::cout << states[i] << "\t"<< k_det << std::endl;
        }    
        k_det *= k_ex_max*phi_f*phi_det;
        arrivals = gsl_ran_poisson (rng, k_det*dt);    
        //std::cout << k_det << "\t" << arrivals <<  std::endl;
        return arrivals;
    }
};

#endif
