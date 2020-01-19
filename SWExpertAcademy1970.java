import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int T;
		T=sc.nextInt();

		for(int test_case = 1; test_case <= T; test_case++)
		{
            int[] money = {50000, 10000, 5000, 1000, 500, 100, 50, 10};
            int x = sc.nextInt();
			System.out.println("#" + test_case);
            for(int i=0; i<money.length;i++)
            {
            	System.out.print(x/money[i] + " ");
                x %= money[i];
            }
            System.out.println();
			

		}
	}
}
