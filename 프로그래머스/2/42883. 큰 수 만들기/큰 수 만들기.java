import java.util.*;
class Solution {
    public String solution(String number, int k) {
        Stack<Integer> stack = new Stack<>();
        for(int i = 0; i < number.length(); i++){
            while(!stack.isEmpty() && number.charAt(stack.peek()) < number.charAt(i) && k > 0){
                stack.pop();
                k--;
            }
            stack.push(i);
        }
        String answer = "";
        while(!stack.isEmpty()){
            Integer idx = stack.pop();
            if(k > 0){
                k--;
                continue;
            }
            answer = number.charAt(idx) + answer;
        }
        return answer;
    }
}