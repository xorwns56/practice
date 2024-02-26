import java.util.*;
class Solution {
    public int solution(String s) {
        int answer = 0;
        char[] chars = s.toCharArray();
        Stack<Character> stack = new Stack<>();
        for(int i = 0; i < chars.length; i++){
            stack.clear();
            for(int j = 0; j < chars.length; j++){
                char c = chars[(i + j) % chars.length];
                if(!stack.isEmpty() && ((stack.peek() == '[' && c == ']') || (stack.peek() == '(' && c == ')') || (stack.peek() == '{' && c == '}'))) stack.pop();
                else stack.push(c);
            }
            if(stack.size() == 0) answer++;
        }
        return answer;
    }
}