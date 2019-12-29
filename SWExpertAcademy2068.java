import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int T = sc.nextInt();
		for(int i=0; i<T; i++)
		{
            int max = sc.nextInt();
            for(int x=0;x<9;x++){
                int input = sc.nextInt();
                if(max<input) max = input;
            }
            System.out.println("#"+(i+1)+" "+max);
		}
	}
}
