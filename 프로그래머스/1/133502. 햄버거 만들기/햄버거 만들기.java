import java.util.*;
class Solution {
    public int solution(int[] ingredient) {
        Stack<Integer> stack = new Stack<>();
        int[] burger_arr = new int[]{1, 2, 3, 1};
        int answer = 0;
        for(int i = 0; i < ingredient.length; i++){
            stack.push(ingredient[i]);
            if(ingredient[i] == burger_arr[burger_arr.length - 1]){
                if(isBurger(stack, burger_arr)) answer++;
            }
        }
        return answer;
    }
    
    public boolean isBurger(Stack stack, int[] burger_arr){
        if(stack.size() < burger_arr.length) return false;
        int idx = burger_arr.length - 1;
        while(idx >= 0){
            if(burger_arr[idx] != (int)stack.peek()) break;
            stack.pop();
            idx--;
        }
        if(idx < 0) return true;
        while(++idx < burger_arr.length) stack.push(burger_arr[idx]);
        return false;
    }
}