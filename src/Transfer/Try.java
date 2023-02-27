package Transfer;

public class Try {
    public static void main(String[] args) {
        Transfer.convertEquatorial(160,0);
        Transfer.convertEquatorial(180,0);

        double[][] a = new double[2][2];
        a[0][0]=150;a[0][1]=0;a[1][0]=120;a[1][1]=0;
        Transfer.convertEquatorial(a);
    }
}
