import java.math.*;

import static java.lang.StrictMath.*;


public class Transfer {

    public static double deg2HMS(double raj ,double decj ,boolean rou){
      //"  convert deg to ra's HMS or dec's DHS'
        String RA,DEC,rs,ds ="";
return 0;
    }

    public static double[] convertEquatorial(float el,float eb){
//''convert Ecliptic(el eb) to Equatorial raj decj'''
        double M_PI =Math.PI;
        double ce = 0.91748213149438;  // Cos epsilon
        double se = 0.39777699580108;  // Sine epsilon
        double dr = M_PI / 180.0;

        double sdec = sin(eb * dr) * ce + cos(eb * dr) * se * sin(el * dr);
        double dec = asin(sdec);  // in radians
        double cos_ra = (cos(el * dr) * cos(eb * dr) / cos(dec));
        double sin_ra = (sin(dec) * ce - sin(eb * dr)) / (cos(dec) * se);
        double ra = atan2(sin_ra, cos_ra);  // in radians

        //# Get RA into range 0 to 360
        if (ra < 0.0) {
            ra = 2 * M_PI + ra;
        }
        else if (ra >2 * M_PI){
            ra = ra - 2 * M_PI;
        }

        double raj = ra*180/M_PI;  //# convert to degrees
        double decj = dec*180/M_PI ;// # convert to degree

        ra = deg2HMS(raj,0,true);
        dec = deg2HMS(0,decj,true);
        double[] NB=new double[2];
        NB[0] = ra;
        NB[1] = dec;

        return NB;
    }



    public static double convertEcliptic(double raj,double decj){
    return 0;
    }
//            '''convert raj deci to elog elat'''
}
