import java.util.Scanner;

public class Solution{
    public static void main(String[] args){
        Scanner input = new Scanner(System.in); //입력을 받기 위한 스캐너
        int t_case = input.nextInt(); //테스트 케이스의 개수 t_case를 받습니다
        for(int i=1;i<=t_case;i++){ //반복문을 t_case만큼 돌립니다
            int month1 = input.nextInt(); //첫 번째 날짜의 월
            int day1 = input.nextInt(); //첫 번째 날짜의 일
            int month2 = input.nextInt(); //두 번째 날짜의 월
            int day2 = input.nextInt(); //두 번째 날짜의 일
            System.out.print("#"+ i + " "); //출력의 각 줄은 # + 테스트 케이스의 번호로 시작합니다
            int sum = month1==month2?day2:maxDay(month1) - (day1 - 1);
            while(++month1<=month2) sum += month1==month2?day2:maxDay(month1);
            System.out.println(sum);
        }
    }
    
   public static int maxDay(int month){ //달의 말일을 반환하는 함수
       switch(month){
           case 4:
           case 6:
           case 9:
           case 11:
               return 30;
           case 2:
               return 28;
           default:
               return 31;
       }
   }
}
