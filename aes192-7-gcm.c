#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define Nb 4
#define BLOCKLEN 16
#define Nk 6
#define Nr 12
#define keyExpSize 208


static const uint8_t sbox[256] = {
  //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };

static const uint8_t rsbox[256] = {
  0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb,
  0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb,
  0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e,
  0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25,
  0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92,
  0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84,
  0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06,
  0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b,
  0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73,
  0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e,
  0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b,
  0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4,
  0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f,
  0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef,
  0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61,
  0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d };

static const uint8_t Rcon[11] = {0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36};
static uint8_t RoundKey[keyExpSize];

uint8_t mul2(uint8_t x){
    uint8_t out = ((x<<1) ^ (((x>>7) & 1) * 0x1b));
    return out;
}

uint8_t mul3(uint8_t x){
    uint8_t out = mul2(x) ^ x;
    return out;
}

uint8_t mul4(uint8_t x){
    uint8_t out = mul2(mul2(x));
    return out;
}

uint8_t mul8(uint8_t x){
    uint8_t out = mul2(mul2(mul2(x)));
    return out;
}

uint8_t mul9(uint8_t x){
    uint8_t out = mul8(x) ^ x;
    return out;
}

uint8_t mulB(uint8_t x){
    uint8_t out = mul8(x) ^ mul2(x) ^ x;
    return out;
}

uint8_t mulD(uint8_t md_in){
    uint8_t out = mul8(md_in) ^ mul4(md_in) ^ md_in;
    return out;
}

uint8_t mulE(uint8_t me_in){
    uint8_t out = mul8(me_in) ^ mul4(me_in) ^ mul2(me_in);
    return out;
}


void State2Matrix(uint8_t state[4][4], const uint8_t input[16]) { 
    int i,j;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            state[j][i] = *input++;  
        }
    }
}

void Matrix2State(uint8_t output[16], uint8_t state[4][4]) { 
    int i,j;
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            *output++ = state[j][i];
        }
    }
}

void key_exp_192(uint8_t key[24]) {
    uint8_t tempa[4];
    for(int i = 0; i < Nk; i++) {

    }

}

void subBytes(uint8_t state[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            state[i][j] = sbox[state[i][j]];
        }
    }
}
void invSubBytes(uint8_t state[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            state[i][j] = rsbox[state[i][j]];
        }
    }
}
void ShiftRows(uint8_t state[4][4]){ 
    int i,j;
    uint8_t temp[4][4];
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
           temp[i][j]=state[i][j];
        }
    }
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
            state[i][j]=temp[i][(j+i)%4]; 
        }
    }
}
void InShiftRows(uint8_t state[4][4]){  
    int i,j;
    uint8_t temp[4][4];
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
           temp[i][j]=state[i][j];
        }
    }
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
            state[i][j]=temp[i][(4+j-i)%4]; 
        }
    }
}

uint8_t XTIME(uint8_t x) {  
	return ((x << 1) ^ ((x & 0x80) ? 0x1b : 0x00));  
}

uint8_t multiply(uint8_t a, uint8_t b) {
	unsigned char temp[8] = { a };
    uint8_t tempmultiply = 0x00;
	int i;
	for (i = 1; i < 8; i++) {
		temp[i] = XTIME(temp[i - 1]); 
	}
	tempmultiply = (b & 0x01) * a;
	for (i = 1; i <= 7; i++) {
		tempmultiply ^= (((b >> i) & 0x01) * temp[i]); 
	}
	return tempmultiply;
}

void MixColumns(uint8_t state[4][4]){
    int i,j;
    uint8_t temp[4][4];
    uint8_t M[4][4]={
        {0x02, 0x03, 0x01, 0x01},
        {0x01, 0x02, 0x03, 0x01},
        {0x01, 0x01, 0x02, 0x03},
        {0x03, 0x01, 0x01, 0x02}
    };
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j){
            temp[i][j] = state[i][j];
        }
    }

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) { 
            state[i][j] = multiply(M[i][0], temp[0][j]) ^ multiply(M[i][1], temp[1][j])
                        ^ multiply(M[i][2], temp[2][j]) ^ multiply(M[i][3], temp[3][j]);
        }
    }
}

void InMixColumns (uint8_t state[4][4]){
    int i,j;
    uint8_t temp[4][4];
    uint8_t M[4][4]={
        {0x0E, 0x0B, 0x0D, 0x09},
        {0x09, 0x0E, 0x0B, 0x0D},
        {0x0D, 0x09, 0x0E, 0x0B},
        {0x0B, 0x0D, 0x09, 0x0E}};  
    for(i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j){
            temp[i][j] = state[i][j];
        }
    }

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            state[i][j] = multiply(M[i][0], temp[0][j]) ^ multiply(M[i][1], temp[1][j])
                        ^ multiply(M[i][2], temp[2][j]) ^ multiply(M[i][3], temp[3][j]);
        }
    }
}

void AddRoundKey(uint8_t state[4][4],uint8_t W[4][4]){  
    int i,j;
    for(i=0;i<4;i++){
        for(j=0;j<4;j++){
            state[j][i]^=W[j][i];
        }
    }
}

static void KeyExpansion16(const uint8_t *Key)
{
  uint32_t i, k;
  uint8_t tempa[4]; 
  
  for (i = 0; i < Nk; ++i)
  {
    RoundKey[(i * 4) + 0] = Key[(i * 4) + 0];
    RoundKey[(i * 4) + 1] = Key[(i * 4) + 1];
    RoundKey[(i * 4) + 2] = Key[(i * 4) + 2];
    RoundKey[(i * 4) + 3] = Key[(i * 4) + 3];
  }

  for (; i < Nb * (Nr + 1); ++i)
  {
    {
      tempa[0]=RoundKey[(i-1) * 4 + 0];
      tempa[1]=RoundKey[(i-1) * 4 + 1];
      tempa[2]=RoundKey[(i-1) * 4 + 2];
      tempa[3]=RoundKey[(i-1) * 4 + 3];
    }

    if (i % Nk == 0)
    {
      {
        k = tempa[0];
        tempa[0] = tempa[1];
        tempa[1] = tempa[2];
        tempa[2] = tempa[3];
        tempa[3] = k;
      }
      {
        tempa[0] = sbox[tempa[0]];
        tempa[1] = sbox[tempa[1]];
        tempa[2] = sbox[tempa[2]];
        tempa[3] = sbox[tempa[3]];
      }

      tempa[0] =  tempa[0] ^ Rcon[i/Nk];
    }
    RoundKey[i * 4 + 0] = RoundKey[(i - Nk) * 4 + 0] ^ tempa[0];
    RoundKey[i * 4 + 1] = RoundKey[(i - Nk) * 4 + 1] ^ tempa[1];
    RoundKey[i * 4 + 2] = RoundKey[(i - Nk) * 4 + 2] ^ tempa[2];
    RoundKey[i * 4 + 3] = RoundKey[(i - Nk) * 4 + 3] ^ tempa[3];
  }
}

static void KeyExpansion(const uint8_t *key,uint8_t w[][4][4]){
    int i,j,r;
    KeyExpansion16(key);
    for(r = 0; r <= Nr; r++) {
        for(i=0;i<4;i++){
            for(j=0;j<4;j++){
                w[r][i][j]=RoundKey[r*16+i+j*4];
            }
        }
    } 
}


void subByte16(uint8_t *RoundText){
    for(int i=0;i<16;i++)
        RoundText[i]=sbox[RoundText[i]];
}
void InvSubByte16(uint8_t *RoundText){
    for(int i=0;i<16;i++)
        RoundText[i]=rsbox[RoundText[i]];
}
void ShiftRow16(uint8_t *RoundText) {
    uint8_t t;
    t=RoundText[1]; RoundText[1]=RoundText[5]; RoundText[5]=RoundText[9]; RoundText[9]=RoundText[13]; RoundText[13]=t;
    t=RoundText[2];RoundText[2]=RoundText[10];RoundText[10]=t;
    t=RoundText[6];RoundText[6]=RoundText[14];RoundText[14]=t;
    t=RoundText[15];RoundText[15]=RoundText[11];RoundText[11]=RoundText[7]; RoundText[7]=RoundText[3]; RoundText[3]=t;
}

void InvShiftRow16(uint8_t *RoundText){
    uint8_t t;
    t=RoundText[13]; RoundText[13]=RoundText[9];RoundText[9]=RoundText[5]; RoundText[5]=RoundText[1]; RoundText[1]=t;
    t=RoundText[2],RoundText[2]=RoundText[10];RoundText[10]=t;
    t=RoundText[6];RoundText[6]=RoundText[14];RoundText[14]=t;
    t=RoundText[3]; RoundText[3]=RoundText[7];RoundText[7]=RoundText[11]; RoundText[11]=RoundText[15]; RoundText[15]=t;
}
void MixColumn16(uint8_t* RoundText){
    uint8_t temp[4];
    for(int i=0;i<4;i++){
        temp[0]=mul2(RoundText[4*i])^mul3(RoundText[1+4*i])^RoundText[2+4*i]^RoundText[3+4*i];
        temp[1]=RoundText[4*i]^mul2(RoundText[1+4*i])^mul3(RoundText[2+4*i])^RoundText[3+4*i];
        temp[2]=RoundText[4*i]^RoundText[1+4*i]^mul2(RoundText[2+4*i])^mul3(RoundText[3+4*i]);
        temp[3]=mul3(RoundText[4*i])^RoundText[1+4*i]^RoundText[2+4*i]^mul2(RoundText[3+4*i]);
        RoundText[4*i]=temp[0];
        RoundText[1+4*i]=temp[1];
        RoundText[2+4*i]=temp[2];
        RoundText[3+4*i]=temp[3];
    }
}

void InvMixColumn16(uint8_t* RoundText) {
    uint8_t temp[4];
    for(int i=0;i<4;i++){
        temp[0]=mulE(RoundText[4*i])^mulB(RoundText[1+4*i])^mulD(RoundText[2+4*i])^mul9(RoundText[3+4*i]);
        temp[1]=mul9(RoundText[4*i])^mulE(RoundText[1+4*i])^mulB(RoundText[2+4*i])^mulD(RoundText[3+4*i]);
        temp[2]=mulD(RoundText[4*i])^mul9(RoundText[1+4*i])^mulE(RoundText[2+4*i])^mulB(RoundText[3+4*i]);
        temp[3]=mulB(RoundText[4*i])^mulD(RoundText[1+4*i])^mul9(RoundText[2+4*i])^mulE(RoundText[3+4*i]);
        RoundText[4*i]=temp[0];
        RoundText[1+4*i]=temp[1];
        RoundText[2+4*i]=temp[2];
        RoundText[3+4*i]=temp[3];
    }
}
void AddRoundKey16(uint8_t* RoundText,int round){
    for(int i=0;i<16;i++)
        RoundText[i]^=RoundKey[round*16+i];
}
void AES192encrypt16(uint8_t* RoundText, uint8_t *key, int R){
    KeyExpansion16(key);
    int round=0;
    AddRoundKey16(RoundText,round);
    for(round=1;round<R;round++){
        subByte16(RoundText);
        ShiftRow16(RoundText);
        MixColumn16(RoundText);
        AddRoundKey16(RoundText,round);
    }
    subByte16(RoundText);
    ShiftRow16(RoundText);
    AddRoundKey16(RoundText,round);
}
void aes_192_enc(uint8_t output[16], const uint8_t key[24], int round, const uint8_t input[16]) {
    int i,j,k;
    uint8_t state[4][4];
    uint8_t w[13][4][4];
    KeyExpansion(key,w);
    State2Matrix(state,input);  
    AddRoundKey(state,w[0]);  
    for(i=1;i<=round;i++){
        subBytes(state);
        ShiftRows(state);
        if(i!=round)
            MixColumns(state);
        AddRoundKey(state,w[i]);
    }
    Matrix2State(output,state);
}

void test_aes_192(){
    //uint8_t in[16]  = {0x6b, 0xc1, 0xbe, 0xe2, 0x2e, 0x40, 0x9f, 0x96, 0xe9, 0x3d, 0x7e, 0x11, 0x73, 0x93, 0x17, 0x2a};
    uint8_t in[16]  = {0xae,0x2d,0x8a,0x57,0x1e,0x03,0xac,0x9c,0x9e,0xb7,0x6f,0xac,0x45,0xaf,0x8e,0x51};
    //8e73b0f7da0e6452c810f32b809079e562f8ead2522c6b7b
    //974104846d0ad3ad7734ecb3ecee4eef
    uint8_t key[24] = { 0x8e, 0x73, 0xb0, 0xf7, 0xda, 0x0e, 0x64, 0x52, 0xc8, 0x10, 0xf3, 0x2b, 0x80, 0x90, 0x79, 0xe5,
                    0x62, 0xf8, 0xea, 0xd2, 0x52, 0x2c, 0x6b, 0x7b};
    const uint8_t out[16] = { 0xbd, 0x33, 0x4f, 0x1d, 0x6e, 0x45, 0xf2, 0x5f, 0xf7, 0x12, 0xa2, 0x14, 0x57, 0x1f, 0xa5, 0xcc};
    uint8_t state[16] = {0xae,0x2d,0x8a,0x57,0x1e,0x03,0xac,0x9c,0x9e,0xb7,0x6f,0xac,0x45,0xaf,0x8e,0x51};
    //aes_192_enc(state, key, 12, in);
    AES192encrypt16(state, key, Nr);
    printf("cipher\n");
    for(int i = 0; i < 16; i++) {
        printf("%02x, ", state[i]);
    }printf("\n");
    
}

void galois_multiply(uint8_t X[16], uint8_t Y[16], uint8_t out[16])
{
    uint8_t V[16];
    uint8_t Z[16];

    for(int i=0;i<16;i++) 
    {
        Z[i]=0x00;
        V[i]=Y[i];
    }

    for (int i = 0; i < 16; i++) 
    {
		for (int j = 0; j < 8; j++) 
		{
			if (X[i] & 1<<(7 - j)) 
			{
				for(int k=0;k<16;k++)
					Z[k]=Z[k]^V[k];
			} 

			if (V[15] & 0x01) 
			{
				for(int k=0;k<16;k++)
				{ 
					if (k!=0)
						if (V[15-k] & 0x01)
							V[16-k] = V[16-k] | 0x80;
					V[15-k]=V[15-k]>>1;
				 }
			  	 V[0]=V[0]^0xe1;
			} 
			else 
			{
				for(int k=0;k<16;k++)
				{ 
					if (k!=0)
						if (V[15-k] & 0x01)
							V[16-k] = V[16-k] | 0x80;
					V[15-k]=V[15-k]>>1;
			       }	
			}
		}
	}
    for(int i=0;i<16;i++)
	{
	  out[i] = Z[i]; 
	}
}

void g_ctr(uint8_t* X, const uint8_t k[24], uint8_t ICB[16], long len_X, int round, uint8_t* out)
{
	uint8_t CB[16];
	uint8_t encrypt_CB[16];
	int inc = 0;
	unsigned int i32count = 1;

	if(len_X > 0)
	{
		for(int i=0;i<16;i++)
			CB[i] = ICB[i];
		for(int i=0;i<len_X;i++)
		{
			i32count++;
			if(i32count==0)
				i32count=1;
            CB[12]=i32count>>24;
            CB[13]=i32count>>16;
            CB[14]=i32count>>8;
            CB[15]=i32count;

            aes_192_enc(encrypt_CB, k, round, CB);

			for(int j =0;j<16;j++)
			{
				out[(i*16)+j] = encrypt_CB[j] ^ X[(i*16)+j];
			}
		}
	}
}

void g_hash(uint8_t H[16], uint8_t* X, unsigned int m, uint8_t out[16])
{
	uint8_t Y[16];

	for( int i=0;i<16;i++)
	{
		Y[i] = 0x00;
	}
	for(int i=0; i < m; i++)
	{
        uint8_t temp[16];
        for(int j =0;j<16;j++)
        {
            temp[j] = Y[j] ^ X[(i*16)+j];
        }
        galois_multiply(temp, H, Y);
    }
    for( int i=0;i<16;i++)
    {
		out[i] = Y[i];
    }
}

void aes192gcm(uint8_t *ciphertext, uint8_t tag[16], uint8_t k[24], const uint8_t IV[12], int round, const uint8_t *plaintext, const unsigned long len_p, const uint8_t* add_data, const unsigned long len_ad) 
{
 
	uint8_t J_0[16]; 
	uint8_t X[16*(len_p+len_ad+1)];
    for(int i = 0; i < 16*(len_p+len_ad+1); i++){
        X[i] = 0;
    }

	for(int i = 0;i<16;i++)
	{
		J_0[i] = 0x00;
		if(i<12)
			J_0[i] = IV[i];
	}
	J_0[15]=0x01;
	uint8_t H[16];
	uint8_t zero[16];
	for(int i = 0;i<16;i++)
	{
		zero[i]=0x00;
	}
	uint8_t P[16*len_p];
	for(int i=0; i<16*len_p; i++)
		P[i]= plaintext[i];

    aes_192_enc(H, k, round, zero);
    g_ctr(P,k,J_0,len_p,round, ciphertext);

	for(int i=0; i<16*len_ad; i++)
		X[i]= add_data[i];
	for(int i=16*len_ad; i<16*(len_ad+len_p); i++)
		X[i]= ciphertext[i-(16*len_ad)];
	for(int i=16*(len_ad+len_p);i< (16*(len_ad+len_p))+7;i++)
	{
		X[i] = ((uint64_t)(16*8*len_ad) >> (8*(7-(i-16*(len_ad+len_p))))) & 0xff ;
	}
	X[(16*(len_ad+len_p))+7] = (16*8*len_ad)  & 0xff;
	for(int i=(16*(len_ad+len_p))+8;i< (16*(len_ad+len_p))+15;i++)
	{
		X[i] = ((uint64_t)(16*8*len_p) >> (8*(7-(i-(16*(len_ad+len_p)))+8))) & 0xff;
	}
	X[(16*(len_ad+len_p))+15] = (16*8*len_p)  & 0xff;
	unsigned long len_x =(len_ad+len_p+1);
	uint8_t Y[16*3];

	g_hash(H, X,len_x, Y);

	uint8_t temp[16*3];
	uint8_t en_counter[16];
    aes_192_enc(en_counter, k, round, J_0);

	for(int i=0;i<16;i++)
	{
		tag[i] = en_counter[i]^Y[i];
	}

}

void key_commit_attack_verify(){
    uint8_t K0[24] = {0x5a,0xdf,0x9b,0xfa,0xed,0xf3,0x75,0x58,0x58,0x35,0x23,0x46,0xd7,0x4d,0x63,0x63,0xa8,0x7c,0xa5,0x90,0xe9,0x78,0x24,0x82};
    uint8_t K1[24] = {0x5a,0xd7,0x9b,0xfa,0xed,0xf3,0x80,0x58,0x58,0x35,0x23,0x46,0xd7,0x4d,0x63,0x63,0xa8,0x7c,0xa5,0x90,0xe9,0x78,0xd1,0x82};
    uint8_t IV[12] = {0x34,0xdf,0x96,0xb6,0x3f,0xb1,0x71,0x5d,0x22,0x9e,0xb9,0x6b};
    uint8_t add_data[16]={0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x01};
    uint8_t plaintext0[48] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x4b,0xd8,0x5,0x59,0x8c,0x87,0x8,0x37,0x1a,0x2f,0xe8,0x29,0x6b,0xa1,0x7,0xf1,0xe4,0x66,0x38,0x9e,0x45,0x26,0x75,0x9d,0x3d,0xae,0xcd,0xb8,0x8a,0xec,0xc5,0x37};
    uint8_t plaintext1[48] = {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0xe6,0xf1,0x2c,0x6e,0x71,0x81,0x3a,0xfe,0x4b,0xcb,0x18,0x42,0xa9,0x89,0xbb,0xe5,0x35,0x26,0x63,0x2b,0xa4,0x3e,0x32,0x1c,0x54,0x14,0x9c,0x4,0x5a,0xb5,0x85,0x1a};
    uint8_t ciphertext0[48] = {0};
    uint8_t ciphertext1[48] = {0};
    uint8_t tag0[16] = {0};
    uint8_t tag1[16] = {0};
    aes192gcm(ciphertext0, tag0, K0, IV, 7, plaintext0, 3, add_data, 1);
    aes192gcm(ciphertext1, tag1, K1, IV, 7, plaintext1, 3, add_data, 1);
    printf("IV\n");
    for(int i = 0; i < 12; i++){
        printf("%02x ",IV[i]);
    }printf("\n");

    printf("AD\n");
    for(int i = 0; i < 16; i++){
        printf("%02x ",add_data[i]);
    }printf("\n");

    printf("plaintext0\n");
    for(int i = 0; i < 48; i++){
        printf("%02x ",plaintext0[i]);
    }printf("\n");
    printf("plaintext1\n");
    for(int i = 0; i < 48; i++){
        printf("%02x ",plaintext1[i]);
    }printf("\n");

    printf("K0\n");
    for(int i = 0; i < 24; i++){
        printf("%02x ",K0[i]);
    }printf("\n");
    printf("K1\n");
    for(int i = 0; i < 24; i++){
        printf("%02x ",K1[i]);
    }printf("\n");

    printf("ciphertext0\n");
    for(int i = 0; i < 48; i++){
        printf("%02x ",ciphertext0[i]);
    }printf("\n");
    printf("ciphertext1\n");
    for(int i = 0; i < 48; i++){
        printf("%02x ",ciphertext1[i]);
    }printf("\n");

    printf("tag0\n");
    for(int i=0;i<16;i++)
	{
		printf("%02x ",tag0[i]);
	}printf("\n");
    printf("tag1\n");
    for(int i=0;i<16;i++)
	{
		printf("%02x ",tag1[i]);
	}printf("\n");

    for(int i = 0; i < 48; i++){
        printf("%x ", ciphertext0[i] ^ ciphertext1[i]);
    }printf("\n");

    for(int i=0;i<16;i++)
	{
		printf("%x ", tag0[i] ^ tag1[i]);
	}printf("\n");

}


void main()
{
    printf("begin verify\n");
    key_commit_attack_verify();
    printf("end verify\n");
}