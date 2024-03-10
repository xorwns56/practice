import java.util.*;
class Solution {
    public int[] solution(int[] numbers) {
        int[] answer = new int[numbers.length];
        for(int i = numbers.length - 1; i >= 0; i--){
            int val = -1;
            for(int j = i + 1; j < numbers.length; j++){
                if(numbers[i] < numbers[j]){
                    val = numbers[j];
                    break;
                }else{
                    if(answer[j] == -1){
                        val = -1;
                        break;
                    }else if(numbers[i] < answer[j]){
                        val = answer[j];
                        break;
                    }
                }
            }
            answer[i] = val;
        }
        return answer;
    }
}