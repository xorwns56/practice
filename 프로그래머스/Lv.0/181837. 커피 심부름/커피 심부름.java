class Solution {
    public int solution(String[] order) {
        int answer = 0;
        for(int i = 0; i < order.length; i++){
            if(order[i].indexOf("cafelatte") < 0) answer += 4500;
            else answer += 5000;
        }
        return answer;
    }
}