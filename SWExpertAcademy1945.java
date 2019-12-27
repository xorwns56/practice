import java.util.Scanner;

public class Solution{
    public static void main(String[] args){
        Scanner input = new Scanner(System.in); //입력을 받기 위한 스캐너
        int t_case = input.nextInt(); //테스트 케이스의 개수 t_case를 받습니다
        for(int i=1;i<=t_case;i++){ //반복문을 t_case만큼 돌립니다
            int t_value = input.nextInt(); //테스트 값을 t_value에 받아옵니다
            System.out.print("#" + i); //출력의 각 줄은 # + 테스트 케이스의 번호로 시작합니다
            int p[] = {2,3,5,7,11}; //배열 p에 소수 5개(2,3,5,7,11)를 담습니다
            for(int j=0;j<p.length;j++){ //p의 크기(5)만큼 반복시킵니다
                int count = 0; //몇 번 나눠지는지 셀 변수 count를 0으로 초기화합니다
                while(t_value % p[j] == 0){
                       t_value /= p[j];
                       count++;
                }//나눠질 때까지 반복시키면서, 나눌 때 마다 그 몫을 t_value에 갱신시키고 count를 증가시킵니다
                System.out.print(" "+count); //소수(p[j]) 별로 나누는 것이 끝날 때마다 count 변수를 출력합니다
            }
            System.out.println(); //테스트 케이스가 끝날 때마다 개행시킵니다
        }
    }
}
