import java.util.Scanner;
 
class Solution
{
    public static void main(String args[]) throws Exception
    {
        Scanner sc = new Scanner(System.in);
        int T = sc.nextInt();
        for(int i=0; i<T; i++)
        {
            int num1 = sc.nextInt();
            int num2 = sc.nextInt();
            char c = num1==num2?'=':num1>num2?'>':'<';
            System.out.println("#"+(i+1)+" "+c);
        }
    }
}
