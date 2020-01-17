import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int x = sc.nextInt();
		for(int i = 0; i <= x; i++)
        {
            System.out.print((int)Math.pow(2, i) + " ");
		}
	}
}
