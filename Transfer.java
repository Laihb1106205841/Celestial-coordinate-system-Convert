package Transfer;

import static java.lang.StrictMath.*;

/**
 * 项目名称：Transfer
 * <p>项目作用：对坐标进行转换。可将赤道坐标系，黄道坐标系，银道坐标系，超星系坐标系进行相互转换<a
 * <p>输入单位：经度（角度制） 纬度（角度制）<a
 *
 *<p>项目版本：1.1
 * <p>项目作者：赖海斌
 *
 * <p>                参考文献：
 * <p><a href="https://www.iau.org/public/themes/naming_stars/">Naming Stars IAU 恒星命名</a>
 * <p><a href="https://docs.astropy.org/en/stable/coordinates/index.html#module-astropy.coordinates">astropy关于 Astronomical Coordinate Systems的介绍</a>
 *<p><a href="https://blog.csdn.net/weixin_43990846/article/details/117912880">天文坐标系转换   作者：Persus</a>
 *<p><a href="https://jingyan.baidu.com/article/3a2f7c2e650c5967aed61169.html">三阶行列式计算方法</a>
 * <p><a href="https://blog.csdn.net/casularm/article/details/302811">赤道和银道坐标转换公式  作者：casularm</a>
 * <p><a href="https://blog.sciencenet.cn/blog-117333-258897.html">坐标变换  作者：钱磊</a>
 * <p><a href="https://zh.wikipedia.org/zh-cn/%E5%9C%B0%E5%B9%B3%E5%9D%90%E6%A8%99%E7%B3%BB">Wiki地平坐标系</a>
 *<p>
 * <p>GitHub：<a href="https://github.com/Laihb1106205841/Celestial-coordinate-system-Convert">Celestial coordinate system convert</a>
 *
 * **/
public class Transfer {
    public final double M_PI=StrictMath.PI;
    public final double deg2rad=M_PI/180.0;
    public double Longitude;
    public double Latitude;

    public Transfer(double Longitude, double Latitude, boolean IsRad){
        this.Longitude=Longitude;
        this.Latitude =Latitude;
        if(IsRad){Longitude=Longitude/deg2rad;Latitude=Latitude/deg2rad;}

    }

    public static double deg2HMS(double ra ,double dec ,boolean rou){
      //"  convert deg to ra's HMS or dec's DHS'
    return ra*dec;
    }

    public static double RadToDegrees(double o){
        return StrictMath.toDegrees(o);
    }
    public static double[]RadToDegrees(double[]o){
        for (double i:o){i=RadToDegrees(i);}
        return o;
    }


    /**
     * <p>赤道坐标系，是一种天球坐标系。过天球中心与地球赤道面平行的平面称为天球赤道面，
     * 它与天球相交而成的大圆称为天赤道。赤道面是赤道坐标系的基本平面。
     * <p>天赤道的几何极称为天极，与地球北极相对的天极即北天极，是赤道坐标系的极。
     * 经过天极的任何大圆称为赤经圈或时圈；与天赤道平行的小圆称为赤纬圈。
     * <p>作天球上一点的赤经圈，从天赤道起沿此赤经圈量度至该点的大圆弧长为纬向坐标，称为赤纬。
     * 赤纬从0°到±90°计量，赤道以北为正，以南为负。赤纬的余角称为极距，从北天极起，从0°到180°计量。
     * **/
    public static double[] convertEquatorial(double el,double eb){
        //el,eb in degrees
//''convert Ecliptic(el eb) to Equatorial raj decj'''
        double M_PI =StrictMath.PI;
        double ce = 0.91748213149438;  // Cos epsilon
        double se = 0.39777699580108;  // Sine epsilon
        double dr = M_PI / 180.0;

        double sdec = sin(eb * dr) * ce + cos(eb * dr) * se * sin(el * dr);//sin deta
        double dec = asin(sdec);  // in radians
        double cos_ra = (cos(el * dr) * cos(eb * dr) / cos(dec));//cos a
        double sin_ra = (sin(dec) * ce - sin(eb * dr)) / (cos(dec) * se);//sin a
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
        //ra = deg2HMS(raj,0,true);
        //dec = deg2HMS(0,decj,true);
        double[] NB=new double[2];
        NB[0] = raj;
        NB[1] = decj;

        System.out.printf("输入黄经："+"%.3f",el);
        System.out.printf(" 输入黄纬："+"%.3f",eb);
        System.out.println();
        System.out.printf("输出赤经："+"%.3f",raj);
        System.out.printf(" 输出赤纬："+"%.3f",decj);
        System.out.println();
        System.out.println();
        return NB;
    }
    public double[] ConvertEquatorial(){return convertEquatorial(Longitude,Latitude);}

    public static double[][] convertEquatorial(double[][]a){
        for(double[] b: a){
            b=convertEquatorial(b[0],b[1]);
        }
        return a;
    }

    /**<p>  正如地图一样，天空中也会有属于自己的坐标，
     * 对应地球上的纬线的在天空被称作赤纬，而经线被称为赤经。
     * <p>  赤经这条线的划定方式便是沿着地球的极轴延伸至无限远，
     * 在直线两端点之间用类似经线的球表面的弧线连成，
     * 赤纬便是将地球的纬线扩大到天空中，而赤道扩张的那条线被叫做天赤道。
     * <p>  而天空中的赤经零度不会再次定在格林尼治天文台的坐标上，
     * 所以人们将春分点（白羊座附近，但经过时间流逝，春分点已经变化位置）
     * 定位赤经0°这样我们可以很好地描述星体在天空中的位置。
     **/
    public static double[] convertEcliptic(double raj,double decj){
        //# rai in degrees
        //# decj in degrees
        double raj1=raj;
        double decj1=decj;
        double M_PI = StrictMath.PI;
        double deg2rad = M_PI / 180.0;
        raj=raj*deg2rad;
        decj=decj*deg2rad;
        double epsilon = 23.439292 * deg2rad;
     /*  epsilon = 23.441884*deg2rad;*/

        double sinb = sin(decj) * cos(epsilon) - cos(decj) * sin(epsilon) * sin(raj);
        double beta = asin(sinb);
        double y = sin(raj) * cos(epsilon) + tan(decj) * sin(epsilon);
        double x = cos(raj);

        double lambdap = atan2(y, x);
        double lambdaa;
        if (lambdap < 0)
        {lambdaa = lambdap + 2 * M_PI;}
        else
        {lambdaa = lambdap;}

        double rlambdaa =lambdaa*180/M_PI;
        double rbeta = beta *180/M_PI;

        double[] NB =new double[2];
        NB[0]=rlambdaa;
        NB[1]=rbeta;
        System.out.printf("输入赤经："+"%.3f",raj1);
        System.out.printf(" 输入赤纬："+"%.3f",decj1);
        System.out.println();
        System.out.printf("输出黄经："+"%.3f",rlambdaa);
        System.out.printf(" 输出黄纬："+"%.3f",rbeta);
        System.out.println();
        System.out.println();

        return NB;
    }
//            '''convert raj deci to elog elat'''

    public static double[][] convertEclptic(double[][]a){
        for(double[] b: a){
            b=convertEcliptic(b[0],b[1]);
        }
        return a;
    }
    public double[] convertEclptic(){return convertEcliptic(Longitude,Latitude);}

/**
 * <p>在任何天球坐标系都需要定义赤道和极点。银道坐标系也一样，需要一条垂直于赤道的子午线作为银经的起点。
 * 经由国际会议决定银道坐标系的银纬和银经分别以“b”和“l”标示，银极的银纬（b）是90°（b=+90°或b=−90°）。
 * 银纬～0°的天体，就位在银河系的盘面（亦即银道面）上，也就是在银河坐标的赤道附近。
 * <p>银道面是整个银道坐标系的基本平面，所有银纬与之平行，银经与之垂直；
 * 银河系成员如恒星、暗物质与气体、尘埃等部绝大部分对称分布在银道面的两侧。
 * 太阳系位处在银道面以北112.7±1.8光年（34.56±0.56秒差距）处，
 * 但因为距离银河系中心30,000光年之遥，相对来说还是非常接近银道面的。
 * <p>银道面和天球的赤道面有123°的夹角，银纬（b）以0°至90°角度为单位度量，
 * 北银极银纬是+90°，位置在后发座，靠近牧夫座的大角星附近；
 * 南银极的银纬是-90°，位置在南天的玉夫座。
 * **/
    public static double[] convertGalactic(double raj,double decj){
        //'''convert raj decj to gl gb'''
        //raj,decj in degrees
        double M_PI = StrictMath.PI;
        double deg2rad = M_PI / 180.0;
        //double gpoleRAJ = 192.85 * deg2rad;
        //double gpoleDECJ = 27.116 * deg2rad;
        double[][] rot =new double[3][3];

        double raj1=raj;
        double decj1=decj;

        raj =raj*deg2rad;
        decj=decj*deg2rad;
        /* Note: Galactic coordinates are defined from B1950 system -
        e.g. must transform from J2000.0 equatorial coordinates to IAU 1958 Galactic coords */
        /* Convert to rectangular coordinates */
        double rx = cos(raj) * cos(decj);
        double ry = sin(raj) * cos(decj);
        double rz = sin(decj);

        /* Now rotate the coordinate axes to correct for the effects of precession */
        /* These values contain the conversion between J2000 and B1950 and from B1950 to Galactic */
        rot[0][0] = -0.054875539726;
        rot[0][1] = -0.873437108010;
        rot[0][2] = -0.483834985808;
        rot[1][0] = 0.494109453312;
        rot[1][1] = -0.444829589425;
        rot[1][2] = 0.746982251810;
        rot[2][0] = -0.867666135858;
        rot[2][1] = -0.198076386122;
        rot[2][2] = 0.455983795705;

        double rx2 = rot[0][0] * rx + rot[0][1] * ry + rot[0][2] * rz;
        double ry2 = rot[1][0] * rx + rot[1][1] * ry + rot[1][2] * rz;
        double rz2 = rot[2][0] * rx + rot[2][1] * ry + rot[2][2] * rz;

    /* Convert the rectangular coordinates back to spherical coordinates */

        double gb = asin(rz2);
        double gl = atan2(ry2, rx2);

        if (gl < 0){gl += 2.0 * M_PI;}

        gl = gl*180/M_PI;
        gb = gb*180/M_PI;

        double[] NB =new double[2];
        NB[0] =gl;
        NB[1] =gb;

        System.out.printf("输入赤经："+"%.3f",raj1);
        System.out.printf(" 输入赤纬："+"%.3f",decj1);
        System.out.println();
        System.out.printf("输出银经："+"%.3f",gl);
        System.out.printf(" 输出银纬："+"%.3f",gb);
        System.out.println();
        System.out.println();

        return NB;
    }

    public static double[][] convertGalactic(double[][]a){
        for(double[] b: a){
            b=convertGalactic(b[0],b[1]);
        }
        return a;
    }
    public double[] convertGalactic(){return convertGalactic(Longitude,Latitude);}
    public static double ThreeMulti (double a,double b,double c){return a*b*c;}
    public static double ThreeMulti (double[][] rot){
        return ThreeMulti(rot[0][0],rot[1][1],rot[2][2])+ThreeMulti(rot[0][1],rot[1][2],rot[2][0])+
                ThreeMulti(rot[0][2],rot[1][0],rot[2][1])-ThreeMulti(rot[0][2],rot[1][1],rot[2][0])-
                ThreeMulti(rot[0][1],rot[1][0],rot[2][2]);
    }
/**
 * <p>北银极的赤道坐标：(alphaGP,deltaGP)=(192.85948,27.12825)单位是度,下同。
 *
 * <p>银心方向l=0,b=0。对应的赤道坐标(alpha,delta)=(266.405,-28.936)。
 *
 * <p>北天极的银经 lCP=122.932。
 *
 * <p>如果天体的赤道坐标为(alpha,delta),则其银道坐标(l,b)可由下面公式求得：
 *
 * <p>反向转换公式如下：
 *
 * <p>sin(delta)=sin(deltaGP)sin(b)+cos(deltaGP)cos(b)cos(lCP-l)
 *
 * <p>cos(delta)sin(alpha-alphaGP)=cos(b)sin(lCP-l)
 *
 * <p>cos(delta)cos(aplha-alphaGP)=cos(deltaGP)sin(b)-sin(deltaGP)cos(b)cos(lCP-l)
 * <p>参考公式出处：https://blog.csdn.net/casularm/article/details/302811
 *
 *<p> <a
 *   href="https://blog.csdn.net/casularm/article/details/302811">{@code 参考公式出处}</a>
 *
 * **/
    public static double[] GalacticConvertEquatorial(double gl,double gb){
        double M_PI = StrictMath.PI;
        double deg2rad = M_PI / 180.0;
        double gpoleRAJ = 192.85948 *deg2rad;
        double gpoleDECJ = 27.12825 *deg2rad;
        double ICP = 122.932 *deg2rad;
        double gl1 =gl;
        double gb1 =gb;
        gl =gl*deg2rad;
        gb =gb*deg2rad;

        double sinDelta =sin(gpoleDECJ)*sin(gb)+cos(gpoleDECJ)*cos(gb)*cos(ICP-gl);
        double Delta = asin(sinDelta);

        double sinAdG =cos(gb)*sin(ICP-gl)/cos(Delta);
        double AdG =asin(sinAdG);
        double Alpha =AdG +gpoleRAJ;

        double eb = Delta;
        double el = Alpha;

        if (el < 0){el += 2.0 * M_PI;}
        el = el*180/M_PI;
        eb = eb*180/M_PI;
        double[] NB =new double[2];
        NB[0] =el;
        NB[1] =eb;

        System.out.printf("输入银经："+"%.3f",gl1);
        System.out.printf(" 输入银纬："+"%.3f",gb1);
        System.out.println();
        System.out.printf("输出赤经："+"%.3f",el);
        System.out.printf(" 输出赤纬："+"%.3f",eb);
        System.out.println();
        System.out.println();
        return NB;
    }
    public static double[][] GalacticConvertEquatorial(double[][]a){
        for(double[] b: a){
            b=GalacticConvertEquatorial(b[0],b[1]);
        }
        return a;
    }
    public double[] GalacticConvertEquatorial(){return GalacticConvertEquatorial(Longitude,Latitude);}


    /**
     *<p>超星系坐标系统的北极点(SGB = 90°)位于银河坐标系统
     * (Galactic Coordinate System)的银经l=47.37°，银纬b=+6.32°之处。
     * <p>在赤道坐标系统中(J2000.0)则是赤经RA=18.9h ，赤纬Dec=+15.7°。
     * 坐标零点(SGB=0°, SGL=0°) 位于银经l= 137.37°，银纬b= 0°之处；
     * 对应到J2000.0分点的赤道坐标系统中大约在赤经RA=2.82h，赤纬 Dec=+59.5°
     * **/
    public static double[] convertSupernova (double raj,double decj){

    //'''convert raj decj to sl sb'''
    // raj,decj in degrees
    double M_PI = StrictMath.PI;
    double deg2rad = M_PI / 180.0;
    double[][] rot =new double[3][3];
    raj=raj*deg2rad;
    decj=decj*deg2rad;

    double rx = cos(raj) * cos(decj);
    double ry = sin(raj) * cos(decj);
    double rz = sin(decj);

        //旋转矩阵
    rot[0][0] = 0.375022041659904;
    rot[0][1] = -0.898320159816291;
    rot[0][2] = 0.228865372515963;
    rot[1][0] = 0.341354889828154;
    rot[1][1] = -0.095717033759227;
    rot[1][2] = -0.935048174501584;
    rot[2][0] = 0.861878940141622;
    rot[2][1] = 0.428787989472616;
    rot[2][2] = 0.270749981762484;

    double rx2 = rot[0][0] * rx + rot[0][1] * ry + rot[0][2] * rz;
    double ry2 = rot[1][0] * rx + rot[1][1] * ry + rot[1][2] * rz;
    double rz2 = rot[2][0] * rx + rot[2][1] * ry + rot[2][2] * rz;

    double gb = asin(rz2);
    double gl = atan2(ry2, rx2);

    if (gl < 0){gl += 2.0 * M_PI;}

    gl = gl*180/M_PI;
    gb = gb*180/M_PI;

    double[] NB =new double[2];
    NB[0] =gl;
    NB[1] =gb;

        System.out.printf("输入赤经："+"%.3f",raj);
        System.out.printf(" 输入赤纬："+"%.3f",decj);
        System.out.println();
        System.out.printf("输出超星系经度："+"%.3f",gl);
        System.out.printf(" 输出超星系纬度："+"%.3f",gb);
        System.out.println();
        System.out.println();

        return NB;
    }
    public static double[][] convertSupernova(double[][]a){
        for(double[] b: a){
            b=convertSupernova(b[0],b[1]);
        }
        return a;
    }
    public double[]convertSupernova(){return convertSupernova(Longitude,Latitude);}
}
