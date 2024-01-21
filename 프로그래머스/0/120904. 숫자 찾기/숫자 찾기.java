class Solution {
    public int solution(int num, int k) {
        int answer = String.valueOf(num).indexOf(String.valueOf(k));
        return answer >= 0 ? answer + 1 : answer;
    }
}