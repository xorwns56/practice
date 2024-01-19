class Solution {
    public int solution(int balls, int share) {
        double answer = 1;
        for(int i = 0; i < share; i++) answer = answer * (balls - i) / (share - i);
        return (int)Math.round(answer);
    }
}