import java.util.*;
class Work{
    int progress;
    int speed;
    Work(int progress, int speed){
        this.progress = progress;
        this.speed = speed;
    }
}
class Solution {
    public int[] solution(int[] progresses, int[] speeds) {
        int day = 0;
        Queue<Work> queue = new LinkedList<>();
        List<Integer> list = new ArrayList<>();
        for(int i = 0; i < progresses.length; i++) queue.add(new Work(progresses[i], speeds[i]));
        while(!queue.isEmpty()){
            int release = 0;
            while(!queue.isEmpty() && queue.peek().progress + queue.peek().speed * day >= 100){
                queue.remove();
                release++;
            }
            if(release > 0) list.add(release);
            day++;
        }
        int[] answer = new int[list.size()];
        for(int i = 0; i < list.size(); i++) answer[i] = list.get(i);
        return answer;
    }
}