import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int x;
        	x = sc.nextInt();
		int sum = 0;
        	for(int i=0;i<4;i++){
        		int rest = x%10;
           		x /= 10;
            		sum += rest;
        	}
        	System.out.println(sum);
	}
}
