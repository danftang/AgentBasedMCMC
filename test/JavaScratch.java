import java.util.HashMap;

public class JavaScratch {

    public static void main(String[] args) {
//        HashMap<Integer,Integer> map1;
        Integer i1 = 1234;
        Integer i2 = i1;

        i2 = 2345;

        System.out.println(i1);
        System.out.println(i2);
    }
}
