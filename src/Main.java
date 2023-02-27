import Transfer.Transfer;

public class Main {
    public static void main(String[] args) {

        Transfer.convertEquatorial(160,0);
        Transfer.convertEquatorial(180,0);

//        double[][] a = new double[2][2];
//        a[0][0]=150;a[0][1]=0;a[1][0]=120;a[1][1]=0;
//        Transfer.convertEquatorial(a);

        Transfer.convertEcliptic(161.534,7.819);
        Transfer.convertEcliptic(0,0);
        Transfer.convertGalactic(160,0);
        Transfer.convertEcliptic(0,0);
        Transfer.convertGalactic(0,0);
        Transfer.convertSupernova(0,0);
        Transfer.GalacticConvertEquatorial(248.152,48.382);
        Transfer a1 =new Transfer(1,1,true);
        a1.ConvertEquatorial();
        Transfer.convertEcliptic(135,150);



    }
}
