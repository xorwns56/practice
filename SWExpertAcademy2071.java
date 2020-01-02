import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int T = sc.nextInt();
		for(int i=0; i<T; i++)
		{
		    int sum = 0;
		    for(int x=0;x<10;x++){
			sum += sc.nextInt();
		    }
		    int avg = (int)(sum/10.0 + 0.5);
		    System.out.println("#"+(i+1)+" "+avg);
		}
	}
}
