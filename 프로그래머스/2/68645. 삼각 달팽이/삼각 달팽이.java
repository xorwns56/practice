class Solution {
    public int[] solution(int n) {
        int[] answer = new int[n * (n + 1) / 2];
        snail(answer, 1, n, 1);
        return answer;
    }
    public void snail(int[] arr, int top, int bot, int num){
        if(top > bot) return;
        int diff = top * top / 2;
        if(top == bot) arr[diff] = num;
        for(int i = top; i < bot && arr[diff] == 0; i++){
            arr[diff] = num++;
            diff += i;
        }
        while(diff + 1 < bot * (bot + 1) / 2 && arr[diff + 1] == 0) arr[diff++] = num++;
        for(int i = bot; i > top && arr[diff] == 0; i--){
            arr[diff] = num++;
            diff -= i;
        }
        snail(arr, top + 2, bot - 1, num);
    }
}