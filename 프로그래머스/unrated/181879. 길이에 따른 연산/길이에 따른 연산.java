class Solution {
    public int solution(int[] num_list) {
        int answer = num_list.length <= 10 ? 1 : 0;
        for(int i = 0; i < num_list.length; i++) answer = num_list.length <= 10 ? (answer * num_list[i]) : (answer + num_list[i]);
        return answer;
    }
}