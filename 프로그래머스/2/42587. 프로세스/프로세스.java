import java.util.*;

class Process{
    int priority;
    boolean answer;
    public Process(int priority, boolean answer){
        this.priority = priority;
        this.answer = answer;
    }
}
class Solution {
    public int solution(int[] priorities, int location) {
        int answer = 0;
        List<Process> queue = new LinkedList<>();
        for(int i = 0; i < priorities.length; i++) queue.add(new Process(priorities[i], i == location));
        while(queue.size() > 0){
            Process curr = queue.remove(0);
            boolean execute = true;
            for(int j = 0; j < queue.size(); j++){
                if(curr.priority < queue.get(j).priority){
                    execute = false;
                    break;
                }
            }
            if(execute){
                answer++;
                if(curr.answer) break;
            }else queue.add(curr);
        }
        return answer;
    }
}