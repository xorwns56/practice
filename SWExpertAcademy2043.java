import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int a = sc.nextInt();
        	int b = sc.nextInt();
		if(a>b) System.out.print(a-b+1);
        	else System.out.print(b-a+1);
	}
}
