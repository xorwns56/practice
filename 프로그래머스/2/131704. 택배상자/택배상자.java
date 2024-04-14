import java.util.*;
class Solution {
    public int solution(int[] order) {
        Stack<Integer> sub_belt = new Stack<>();
        Queue<Integer> main_belt = new LinkedList<>();
        for(int i = 0; i < order.length; i++) main_belt.add(i + 1);
        int order_idx = 0;
        while(!main_belt.isEmpty()){
            int box = main_belt.remove();
            sub_belt.push(box);
            while(order_idx < order.length && !sub_belt.isEmpty() && order[order_idx] == sub_belt.peek()){
                sub_belt.pop();
                order_idx++;
            }
        }
        return order_idx;
    }
}