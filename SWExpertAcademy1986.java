import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int T = sc.nextInt();
		int[] sum = new int[T];
		for(int i=0; i<T; i++)
		{
            int n = sc.nextInt();
            sum[i]=0;
        	for(int j=1;j<=n;j++)	sum[i]+=j%2==1?j:-j;
            
            System.out.println("#"+(i+1)+" "+sum[i]);
		}
	}
}
