import java.util.Scanner;
import java.util.Arrays;

class Solution
{
	public static void main(String args[]) throws Exception
	{
    		Scanner sc = new Scanner(System.in);
		int N = sc.nextInt();
		int array[] = new int[N];
    		for(int i=0;i<N;i++) array[i] = sc.nextInt();
   		Arrays.sort(array);
    		System.out.println(array[N/2]);
	}
}
