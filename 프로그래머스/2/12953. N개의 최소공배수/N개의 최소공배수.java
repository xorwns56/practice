import java.util.*;
class Solution {
    public int solution(int[] arr) {
        Stack<Integer> stack = new Stack<>();
        for(int i = 0; i < arr.length; i++){
            if(stack.isEmpty()) stack.push(arr[i]);
            else stack.push(lcm(stack.pop(), arr[i]));
        }
        return stack.pop();
    }
    public int gcd(int a, int b){
        if(b == 0) return a;
        return gcd(b, a % b);
    }
    public int lcm(int a, int b){
        return a * b / gcd(a, b);
    }
}