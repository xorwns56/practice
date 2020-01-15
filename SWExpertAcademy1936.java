import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int a = sc.nextInt();
		int b = sc.nextInt();
		if(a-b<0){
		    if(b-a == 1) System.out.println("B");
		    else System.out.println("A");
		}else{
		    if(a-b == 1)System.out.println("A");
		    else System.out.println("B");
		}
        
	}
}
