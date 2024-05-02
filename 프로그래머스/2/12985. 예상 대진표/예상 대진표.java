class Solution
{
    public int solution(int n, int a, int b)
    {
        int[] player = new int[n];
        for(int i = 0; i < n; i++) player[i] = i + 1;
        return game_ab(player, a, b, 1);
    }
    public int game_ab(int[] player, int a, int b, int round){
        int[] winner = new int[player.length / 2];
        for(int i = 0; i < player.length / 2; i++){         
            winner[i] = player[i * 2];
            if((player[i * 2] == a && player[i * 2 + 1] == b) || (player[i * 2] == b && player[i * 2 + 1] == a)) return round;
            else if(player[i * 2 + 1] == a || player[i * 2 + 1] == b) winner[i] = player[i * 2 + 1];
        }
        return game_ab(winner, a, b, round + 1);
    }
}