import java.util.*;
class Solution {
    public int solution(int[] order) {
        Stack<Integer> sub_belt = new Stack<>();
        int box = 1;
        int order_idx = 0;
        while(box <= order.length){
            sub_belt.push(box++);
            while(order_idx < order.length && !sub_belt.isEmpty() && order[order_idx] == sub_belt.peek()){
                sub_belt.pop();
                order_idx++;
            }
        }
        return order_idx;
    }
}