class Solution {
    public String solution(int[] numbers, String hand) {
        String answer = "";
        int l = 9;
        int r = 11;
        for(int i = 0; i < numbers.length; i++){
            int pos = numbers[i] != 0 ? numbers[i] - 1 : 10;
            int x = pos / 3;
            int y = pos % 3;
            boolean left = false;
            if(y == 0) left = true;
            else if(y == 2) left = false;
            else{
                int l_diff = Math.abs(l / 3 - x) + Math.abs(l % 3 - y);
                int r_diff = Math.abs(r / 3 - x) + Math.abs(r % 3 - y);
                if(l_diff == r_diff) left = hand.equals("left");
                else if(l_diff < r_diff) left = true;
                else left = false;
            }
            if(left){
                answer += "L";
                l = pos;
            }else{
                answer += "R";
                r = pos;
            }
        }
        return answer;
    }
}