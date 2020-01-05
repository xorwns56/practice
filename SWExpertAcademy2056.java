import java.util.Scanner;

class Solution
{
	public static void main(String args[]) throws Exception
	{
		Scanner sc = new Scanner(System.in);
		int T = sc.nextInt();
		for(int i=0; i<T; i++)
		{
		    String input = sc.next();
		    String year = input.substring(0, 4);
			String month = input.substring(4, 6);
			String day = input.substring(6, 8);
			
			System.out.print("#"+(i+1)+" ");
		    if(isCorrectMonth(month)&&isCorrectDay(month,day)) System.out.println(year + "/" + month + "/" + day);
		    else System.out.println(-1);
		}
	}
    
    private static boolean isCorrectMonth(String month){
        int m = Integer.parseInt(month);
        return 1 <= m && m <= 12;
    }
    
    private static boolean isCorrectDay(String month, String day){
        int m = Integer.parseInt(month);
        int d = Integer.parseInt(day);
        switch(m){
            case 1:
            case 3:
            case 5:
            case 7:
            case 8:
            case 10:
            case 12:
                return d<=31;
            case 4:
            case 6:
            case 9:
            case 11:
                return d<=30;
            case 2:
                return d<=28;
            default:
                return false;
        }
    }
}
