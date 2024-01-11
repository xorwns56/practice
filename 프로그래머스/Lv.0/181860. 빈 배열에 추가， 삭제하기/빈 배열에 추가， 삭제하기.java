import java.util.*;
class Solution {
    public Stack<Integer> solution(int[] arr, boolean[] flag) {
        Stack<Integer> stk = new Stack<>();
        for(int i = 0; i < arr.length; i++){
            if(flag[i]){
                for(int j = 0; j < arr[i] * 2; j++) stk.push(arr[i]);
            }else{
                for(int j = 0; j < arr[i]; j++) stk.pop();
            }
        }
        return stk;
    }
}