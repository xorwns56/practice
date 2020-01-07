import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		String s=sc.next();
		char[] c = s.toCharArray();
        for(int i=0;i<c.length;i++){
            if('a'<=c[i]&&c[i]<='z') c[i] = (char)(c[i] - 'a' + 'A');
            System.out.print(c[i]);
        }
	}
}
