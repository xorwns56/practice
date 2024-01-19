class Solution {
    public int[] solution(int[] numbers, String direction) {
        int[] answer = new int[numbers.length];
        int shift = direction.equals("left") ? 1 : -1;
        for(int i = 0; i < answer.length; i++){
            int new_i = i + shift;
            if(new_i < 0) new_i += numbers.length;
            else if(new_i >= numbers.length) new_i %= numbers.length;
            answer[i] = numbers[new_i];
        }
        return answer;
    }
}