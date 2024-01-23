class Solution {
    public String[] solution(String my_str, int n) {
        String[] answer = new String[(my_str.length() + n - 1) / n];
        for(int i = 0; i < answer.length; i++) answer[i] = my_str.substring(n * i, Math.min(my_str.length(), n * (i + 1)));
        return answer;
    }
}