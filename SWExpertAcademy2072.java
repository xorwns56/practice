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
			    int input = sc.nextInt();
			    if((input&1)==1) sum += input;
		    }
		    System.out.println("#"+(i+1)+" "+sum);
		}
	}
}
