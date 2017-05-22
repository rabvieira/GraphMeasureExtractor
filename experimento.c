/*
-----------------------------------------------------------------------------
 Nome      : experimento.c
 Autor     : Renan A. B. Vieira
 Descriçao : Estudo de diferentes tecnicas em diversos datasets
 Consideracoes atuais: Aplica-se apenas uma medida de similaridade por vez;
                       Sigma, por ora, esta constante.
-----------------------------------------------------------------------------
                                            DADOS DE SAIDA:
sigma // normalizou {1->sim/0->nao} //distancia utilizada // %pre-rotulados  // e/k// grau medio // %LGC // %FH
*/
    /*bibliotecas*/
#include <igraph.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
    /*constantes*/
#define normaliza 1       //para normalizar, normaliza = 1
#define alfa 0.5
#define sigma_GLOBAL 0.5
#define distance 0        //0 -> Euclidiana / 1 -> Pearson / 2 -> Distancia sugerida no paper de Funcoes Harmonicas
#define exec 100            //execucoes de cada configuracao 100
#define rot_inicio 0.01   // % de pre-rotulados (inicio) 0.01
#define rot_fim    0.2   // % de pre-rotulados (fim)    0.2
//limiar de corte
#define e_inicio 0.0
#define e_fim    1.00
//K vizinhos mais proximos
#define k_inicio 1
#define k_fim    20

/*
//Imprime matriz
void print_matrix(igraph_matrix_t *m, FILE *f, int mode){
    long int i, j;
    for (i=0 ; i<igraph_matrix_nrow(m) ; i++){
        for (j=0 ; j<igraph_matrix_ncol(m) ; j++){
            if(mode)
                fprintf(f, "%li\t", (long int)MATRIX(*m, i, j));
            else
                fprintf(f, "%.3f\t", (float)MATRIX(*m, i, j));
        }
        fprintf(f, "\n");
    }
}

//Imprime vetor
void print_vector(igraph_vector_t *v, FILE *f){
    long int i;
    for (i=0 ; i<igraph_vector_size(v) ; i++)
        fprintf(f, "%.3f\t", (float)VECTOR(*v)[i]);
}

//Armazena a matriz em um arquivo texto
void armazena(igraph_matrix_t *m){
    int i,j;
    FILE *outstream;
    outstream = fopen("matrix X","w");
    for(i=0 ; i<igraph_matrix_nrow(m) ; i++){
        for(j=0 ; j<igraph_matrix_ncol(m) ; j++){
            fprintf(outstream,"%.3f  ",MATRIX(*m,i,j));
        }
        fprintf(outstream,"\n");
    }
    fclose(outstream);
}
*/

//Salva todos os dados em um arquivo texto
void armazena_dados(char title[], igraph_matrix_t *dados){
    int i,j;
    FILE *outstream;
    outstream = fopen(title,"w");
    for(i=0 ; i<igraph_matrix_nrow(dados) ; i++){
        for(j=0 ; j<igraph_matrix_ncol(dados) ; j++){
            fprintf(outstream,"%.2f  ",MATRIX(*dados,i,j));
        }
        fprintf(outstream,"\n");
    }
    fclose(outstream);
}

//Constroi a matriz com os dados: K/e, % de pre-rotulados, % de classificacao obtida com os diversos algoritmos
void armazena_parametros(igraph_matrix_t *dados, float par, float labeled, int lin, float km, float taxa_LGC,float taxa_FH){
    MATRIX(*dados, lin, 0) = sigma_GLOBAL; //sigma utilizado
    MATRIX(*dados, lin, 1) = normaliza;    //informa se os dados foram normalizados ou nao
    MATRIX(*dados, lin, 2) = distance;     //informa a distancia utilizada
    MATRIX(*dados, lin, 3) = labeled *100; //porcentagem de rotulados de cada classe
    MATRIX(*dados, lin, 4) = par;          //parametro e (corte e) ou k (KNN)
    MATRIX(*dados, lin, 5) = km;           //grau medio da rede
    MATRIX(*dados, lin, 6) = taxa_LGC*100; //porcentagem de acertos aplicando LGC
    MATRIX(*dados, lin, 7) = taxa_FH*100;  //porcentagem de acertos aplicando FH
    igraph_matrix_add_rows(dados,1);
}

//Encontra o indice do maior elemento da linha lin
int max_row(igraph_matrix_t *m, int lin){
    int j,index=0;
    float max;
    max = MATRIX(*m,lin,0);
    for(j=0 ; j<igraph_matrix_ncol(m) ; j++){
        if(MATRIX(*m,lin,j)>max){
            max = MATRIX(*m,lin,j);
            index = j;
        }
    }
    return (index);
}

//Multiplicacao de matrizes
int mul_matrix(igraph_matrix_t *A, igraph_matrix_t *B, igraph_matrix_t *res){
    long int i,j,k;
    if(igraph_matrix_ncol(A)!=igraph_matrix_nrow(B))
        return (-1);
    igraph_matrix_init(res, igraph_matrix_nrow(A), igraph_matrix_ncol(B));
    igraph_matrix_null(res);
    for(i=0 ; i<igraph_matrix_nrow(res) ; i++)
        for(j=0 ; j<igraph_matrix_ncol(res) ; j++)
            for(k=0 ; k<igraph_matrix_ncol(A) ; k++)
                MATRIX(*res, i, j) += MATRIX(*A, i, k) * MATRIX(*B, k, j);
    return (0);
}

//Algoritmo numerico para inversao de matriz
void inverse_matrix(igraph_matrix_t *m, igraph_matrix_t *inv){
    float ratio,a;
    int i, j, k, n;
    igraph_matrix_init(inv, igraph_matrix_nrow(m), 2*igraph_matrix_nrow(m));
    igraph_matrix_null(inv);
    n = (int)igraph_matrix_nrow(m);
//copia a matriz m em inv
    for(i=0 ; i<n ; i++)
        for(j=0 ; j<n ; j++)
            MATRIX(*inv,i,j) = MATRIX(*m,i,j);
//Algoritmo Numerico
    for(i=0 ; i<n ; i++){
        for(j=n ; j<2*n ; j++){
            if(i==(j-n))
                MATRIX(*inv,i,j) = 1.0;
            else
                MATRIX(*inv,i,j) = 0.0;
        }
    }
    for(i=0 ; i<n ; i++){
        for(j=0 ; j<n ; j++){
            if(i!=j){
                ratio = MATRIX(*inv,j,i)/MATRIX(*inv,i,i);
                for(k=0 ; k<2*n ; k++){
                    MATRIX(*inv,j,k) -= ratio * MATRIX(*inv,i,k);
                }
            }
        }
    }
    for(i=0 ; i<n ; i++){
        a = MATRIX(*inv,i,i);
        for(j=0; j<2*n; j++){
            MATRIX(*inv,i,j) /= a;
        }
    }
//Transporta a inversa para inv.
    igraph_matrix_null(m);
    for(i=0 ; i<n ; i++)
        for(j=n ; j<2*n ; j++)
            MATRIX(*m,i,j-n) = MATRIX(*inv,i,j);
    for(j=0 ; j<n ; j++)
        igraph_matrix_remove_col(inv,j);
    igraph_matrix_null(inv);
    for(i=0 ; i<n ; i++)
        for(j=0 ; j<n ; j++)
            MATRIX(*inv,i,j) = MATRIX(*m,i,j);
}

//Soma da i-esima linha da matriz W
float sum_W(igraph_matrix_t *W, int i){
    int k;
    float sum=0;
    for(k=0 ; k<igraph_matrix_ncol(W) ; k++){
        sum += MATRIX(*W, i, k);
    }
    return (sum);
}

//Algoritmo "Semi-Supervised Using Gaussian Fields and Harmonic Functions" e retorna a taxa de classificacao
float Fharmonicas(igraph_vector_t *gabarito, igraph_vector_t *rotulos, igraph_matrix_t *m, int classes, int l){
    int i,j,k, lab=0, unl=0, acertos=0;
    igraph_vector_t D;             //vetor D' (soma das linhas da matriz de adj.)
    igraph_vector_t gab_aux;
    igraph_vector_t rot_aux;
    igraph_matrix_t matrix_Yl;     //matriz de rotulos labeled
    igraph_matrix_t matrix_X;      //matriz X = I - Puu (particao UU da matriz de transicao de probabilidades)
    igraph_matrix_t matrix_inv_X;  //matriz inversa de X
    igraph_matrix_t matrix_Pul;    //matriz Pul (particao UL da matriz de transicao de probabilidades)
    igraph_matrix_t matrix_prod;   //produto Pul*Yl
    igraph_matrix_t matrix_Hu;     //matriz Hu = X(inv) * Pul*Yl
    igraph_vector_init(&gab_aux, igraph_vector_size(gabarito));
    igraph_vector_init(&rot_aux, igraph_vector_size(rotulos));
    igraph_vector_init(&D, l);
//copia o gabarito e os rotulos, pois estes sao modificaveis neste  algoritmo
    for(i=0 ; i<igraph_vector_size(gabarito) ; i++){
        VECTOR(gab_aux)[i] = VECTOR(*gabarito)[i];
        VECTOR(rot_aux)[i] = VECTOR(*rotulos)[i];
        VECTOR(D)[i] = sum_W(m,i);
    }
//conta quantos sao os rotulados e os nao-rotulados
    for(i=0 ; i<igraph_vector_size(rotulos) ; i++)
        if((int)VECTOR(*rotulos)[i] != -1)
            lab++;
        else
            unl++;
//Gera a matriz de rotulos (lab x classes) e organiza as particoes
    igraph_matrix_init(&matrix_Yl, lab, classes);
    igraph_matrix_null(&matrix_Yl);
    j=0;
    for(i=0 ; i<igraph_vector_size(rotulos) ; i++){
        if((int)VECTOR(*rotulos)[i] != -1){
            MATRIX(matrix_Yl, j, (int)VECTOR(*rotulos)[i]) = 1;
            igraph_matrix_swap_rows(m,i,j);
            k = VECTOR(rot_aux)[i];
            VECTOR(rot_aux)[i] = VECTOR(rot_aux)[j];
            VECTOR(rot_aux)[j] = k;
            k = VECTOR(gab_aux)[i];
            VECTOR(gab_aux)[i] = VECTOR(gab_aux)[j];
            VECTOR(gab_aux)[j] = k;
            k = VECTOR(D)[i];
            VECTOR(D)[i] = VECTOR(D)[j];
            VECTOR(D)[j] = k;
            j++;
        }
    }
//Gera a matrix X (X = I - Wij/di) para i diferente de j e 1 caso contrario
    igraph_matrix_init(&matrix_X, unl, unl);
    for(i=0 ; i<unl ; i++){
        for(j=0 ; j<unl ; j++){
            if(i==j)
                MATRIX(matrix_X, i, j) = 1;
            else
                MATRIX(matrix_X, i, j) = -MATRIX(*m,i+lab,j+lab) / VECTOR(D)[i+lab];
        }
    }
//Inverte a matrix X
    inverse_matrix(&matrix_X,&matrix_inv_X);
    igraph_matrix_destroy(&matrix_X);
//Monta a matriz Pul
    igraph_matrix_init(&matrix_Pul, unl, lab);
    for(i=0 ; i<unl ; i++)
        for(j=0 ; j<lab ; j++)
            MATRIX(matrix_Pul, i, j) = MATRIX(*m,i+lab,j) / VECTOR(D)[i+lab];
//Multiplica Pul*Yl
    mul_matrix(&matrix_Pul,&matrix_Yl,&matrix_prod);
    igraph_matrix_destroy(&matrix_Pul);
    igraph_matrix_destroy(&matrix_Yl);
//Multiplica X(inv)*Pul*Yl
    mul_matrix(&matrix_inv_X,&matrix_prod,&matrix_Hu);
    igraph_matrix_destroy(&matrix_inv_X);
    igraph_matrix_destroy(&matrix_prod);
//Obtem a taxa de classificacao
    acertos = lab; //dados pre-rotulados
    for (i=lab ; i<igraph_vector_size(gabarito) ; i++)
        if(max_row(&matrix_Hu,i-lab) == (int)VECTOR(gab_aux)[i])
            acertos++;
    igraph_matrix_destroy(&matrix_Hu);
    igraph_vector_destroy(&gab_aux);
    igraph_vector_destroy(&rot_aux);
    igraph_vector_destroy(&D);
    return ((float)acertos/(float)igraph_vector_size(gabarito));
}

//Aplica o algoritmo "Local and Global Consistency" e retorna a taxa de classificacao
float local_global_consistency(igraph_vector_t *gabarito, igraph_vector_t *rotulos, igraph_matrix_t *m, int classes, int l){
    int i,j, acertos=0;
    igraph_matrix_t matrix_Y;       //matriz de rotulos
    igraph_matrix_t matrix_X;       //matriz X = I - alfa*S
    igraph_matrix_t matrix_inv_X;   //matriz inversa de X
    igraph_matrix_t matrix_prod;    //produto X_inv * Y
    igraph_vector_t D;              //vetor D' (soma das linhas da matriz de adj.)
//Gera o vetor D'
    igraph_vector_init(&D, l);
    for(i=0 ; i<l ; i++)
        VECTOR(D)[i] = sum_W(m,i);
//Gera a matriz de rotulos (l x classes)
    igraph_matrix_init(&matrix_Y, l, classes);
    igraph_matrix_null(&matrix_Y);
    for(i=0 ; i<igraph_vector_size(rotulos) ; i++)
        if((int)VECTOR(*rotulos)[i] != -1)
            MATRIX(matrix_Y, i, (int)VECTOR(*rotulos)[i]) = 1;
//Gera a matrix X (X = I - alfa*S), onde {S = Wij / sqrt(Dii)*sqrt(Djj)} para i diferente de j e 1 caso contrario
    igraph_matrix_init(&matrix_X, l, l);
    for(i=0 ; i<l ; i++){
        for(j=0 ; j<l ; j++){
            if(i==j)
                MATRIX(matrix_X, i, j) = 1;
            else
                MATRIX(matrix_X, i, j) = ((-alfa * MATRIX(*m, i, j)) / (sqrt(VECTOR(D)[i]) * sqrt(VECTOR(D)[j])));
        }
    }
//Inverte a matrix X
    inverse_matrix(&matrix_X,&matrix_inv_X);
    igraph_matrix_destroy(&matrix_X);
//Multiplica X_inv * Y
    mul_matrix(&matrix_inv_X,&matrix_Y,&matrix_prod);
    igraph_matrix_destroy(&matrix_inv_X);
    igraph_matrix_destroy(&matrix_Y);
//Obtem a taxa de classificacao
    for (i=0 ; i<igraph_vector_size(gabarito) ; i++)
        if(max_row(&matrix_prod,i) == (int)VECTOR(*gabarito)[i])
            acertos++;
    igraph_matrix_destroy(&matrix_prod);
    igraph_vector_destroy(&D);
    return ((float)acertos/(float)igraph_vector_size(gabarito));
}

//Seleciona labeled % de exemplos (de cada classe) pre-rotulados. Obs: Pega o teto da %
void selecionar_pre_rotulados(igraph_vector_t *gabarito, igraph_vector_t *rotulos, float labeled, int classes){
    int i, j, dado;
    igraph_vector_t cont;
    igraph_vector_init(&cont,classes);
    igraph_vector_init(rotulos,igraph_vector_size(gabarito));
    for(i=0 ; i<igraph_vector_size(gabarito) ; i++)
        VECTOR(*rotulos)[i] = -1;                   //desrotula o vetor de rotulos
    for(i=0 ; i<igraph_vector_size(gabarito) ; i++)
        VECTOR(cont)[(int)VECTOR(*gabarito)[i]] ++; //conta quantos exemplos tem em cada classe
//Determina quantos exemplos pre-rotular em cada classe
    for(i=0 ; i<igraph_vector_size(&cont) ; i++){
        VECTOR(cont)[i] = (labeled*VECTOR(cont)[i])+0.5; //arredonda pra cima
        if(VECTOR(cont)[i]<1)                            //garante que pelo menos 1 exemplo de cada classe
            VECTOR(cont)[i] = 1;
    }
//Pre-rotula os exemplos de forma randomica
    for(i=0 ; i<igraph_vector_size(&cont) ; i++){
        if(VECTOR(cont)[i]!=0){   //garante a convergencia para datasets que nao possuem exemplos em determinada classe
            do{
                j=(int)rand()%igraph_vector_size(gabarito);
                if(VECTOR(*gabarito)[j]==i && VECTOR(*rotulos)[j]==-1){
                    VECTOR(*rotulos)[j] = i;
                    VECTOR(cont)[i] --;
                }
            }while((int)VECTOR(cont)[i]>0);
        }
    }
    igraph_vector_destroy(&cont);
}

//Calcula o grau medio da rede
float grau_medio(igraph_matrix_t *m){
    int i,sum=0;
    for(i=0 ; i<igraph_matrix_nrow(m) ; i++)
        sum += sum_W(m,i);
    return((float)sum/igraph_matrix_nrow(m));
}

//Procura o indice que possui determinado peso w
int index_v(igraph_vector_t *v, float w){
    int i;
    for (i=0 ; i<igraph_vector_size(v) ; i++)
        if((float)VECTOR(*v)[i] == w)
            return(i);
}

//Corta as conexoes abaixo do limiar e
void gerar_rede_e(igraph_matrix_t *W, igraph_matrix_t *m, float e){
    int i,j;
    igraph_matrix_init(m, igraph_matrix_nrow(W), igraph_matrix_ncol(W));
    igraph_matrix_null(m);
    for (i=0 ; i<igraph_matrix_nrow(m) ; i++)
        for (j=0 ; j<igraph_matrix_ncol(m) ; j++)
            if(MATRIX(*W,i,j)>=e)
                MATRIX(*m,i,j) = MATRIX(*W,i,j);
}

//Conecta os K vizinhos mais proximos
void gerar_rede_KNN(igraph_matrix_t *W, igraph_matrix_t *m, int k){
    int i,j;
    igraph_vector_t v,vsort;
    igraph_vector_init(&v,igraph_matrix_ncol(W));
    igraph_vector_init(&vsort,igraph_matrix_ncol(W));
    igraph_matrix_init(m, igraph_matrix_nrow(W), igraph_matrix_ncol(W));
    igraph_matrix_null(m);
    for (i=0 ; i<igraph_matrix_nrow(m) ; i++){
        igraph_matrix_get_row(W,&v,i);
        igraph_matrix_get_row(W,&vsort,i);
        igraph_vector_sort(&vsort); //ordena
        for(j=0 ; j<k ; j++)
            MATRIX(*m,i,index_v(&v,VECTOR(vsort)[igraph_matrix_ncol(W)-1 -j])) = 1;
            //MATRIX(*m,i,index_v(&v,VECTOR(vsort)[igraph_matrix_ncol(W)-1 -j])) = VECTOR(vsort)[igraph_matrix_ncol(W)-1 -j];
    }
//simetriza o grafo -> nao-direcionado
    for(i=0 ; i<igraph_matrix_nrow(m) ; i++)
        for (j=0 ; j<igraph_matrix_ncol(m) ; j++)
            if(MATRIX(*m,i,j)<MATRIX(*m,j,i))
                MATRIX(*m,i,j) = MATRIX(*m,j,i);
    igraph_vector_destroy(&v);
    igraph_vector_destroy(&vsort);
}

//Calcula a media do atributo referente a coluna col
float media(igraph_matrix_t *m, int col){
    int i;
    float res=0;
    for(i=0 ; i<igraph_matrix_nrow(m) ; i++)
        res += MATRIX(*m, i, col);
    return(res/igraph_matrix_nrow(m));
}

//Calcula o desvio padrao do atributo referente a coluna col
float desvio_padrao(igraph_matrix_t *m, int col, float media){
    int i;
    float res=0;
    for(i=0 ; i<igraph_matrix_nrow(m) ; i++)
        res += pow((MATRIX(*m, i, col) - media),2);
    res = res/(igraph_matrix_nrow(m)-1);
    return (sqrt(res));
}

//Calcula a similaridade entre os exemplos
float similarity(igraph_matrix_t *m, int i, int j, int ini, int fim){
    float d=0, sigma=0, media_k=0, desvp_k=0;
    int k;
    sigma = sigma_GLOBAL;
    switch(distance){
        case 0: //Euclidian Distance
            for(k=ini ; k<fim ; k++)
                d += pow(MATRIX(*m, i, k) - MATRIX(*m, j, k),2);
            return exp(-d / (2 * pow(sigma,2)));
        case 1: //Pearson Distance
            for(k=ini ; k<fim ; k++){
                media_k = media(m,k);
                desvp_k = desvio_padrao(m,k,media_k);
                d += pow((((MATRIX(*m, i, k) - media_k)/desvp_k) - ((MATRIX(*m, j, k) - media_k)/desvp_k)),2);
            }
            d = sqrt(d);
            return exp(-d / (2 * pow(sigma,2)));
        case 2: //F. H. Distance
            for(k=ini ; k<fim ; k++){
                media_k = media(m,k);
                desvp_k = desvio_padrao(m,k,media_k);
                d += (pow(MATRIX(*m,i,k) - MATRIX(*m,j,k),2)) / (pow(desvp_k,2));
            }
            return exp(-d / (2 * pow(sigma,2)));
    }
    return(-1); //se modo de similaridade invalido
}

//Normaliza os atributos (colunas) [0,1]
void normaliza_dados(igraph_matrix_t *dataset, int col){
    int i;
    float max, min;
    max = MATRIX(*dataset,0,col);
    min = max; //Apenas para inicializar
    for(i=0 ; i<igraph_matrix_nrow(dataset) ; i++){
        if(MATRIX(*dataset,i,col)>max)
            max = MATRIX(*dataset,i,col);
        if(MATRIX(*dataset,i,col)<min)
            min = MATRIX(*dataset,i,col);
    }
    for(i=0 ; i<igraph_matrix_nrow(dataset) ; i++)
        MATRIX(*dataset,i,col) = (MATRIX(*dataset,i,col) - min) / (max - min);
}

//Extrai do dataset as informacoes relevantes e que sao fixas durante o programa
int gerar_matrizes_fixas(igraph_vector_t *gabarito, igraph_matrix_t *matrix_W, char title[], int l, int c, int rotulo_mode,
                         int separador){
    int ini, fim, i, j;
    float data;
    igraph_matrix_t dataset;  //dataset
    FILE *instream;
    instream = fopen(title,"r");
    if (instream == NULL){
        printf("Erro ao abrir o arquivo para leitura\n");
        return (-1);
    }
//Verifica se os rotulos estao na primeira ou ultima coluna do dataset e configura o range de atributos
    if(rotulo_mode!=0)
        rotulo_mode = c-1;
    if(rotulo_mode==0){
        ini = 1;
        fim = c;
    }
    else{
        ini = 0;
        fim = c-1;
    }
//Armazena o dataset em uma matriz (l x c)
    igraph_matrix_init(&dataset, l, c);
    for(i=0 ; i<l ; i++){
        for(j=0 ; j<c ; j++){
            if(separador) //separado por virgula
                fscanf(instream,"%f,",&data);
            else          //separado por espaco
                fscanf(instream,"%f",&data);
            MATRIX(dataset, i, j) = data;
        }
    }
//Normaliza o dataset
    if(normaliza){
        for(j=ini ; j<fim ; j++)
            normaliza_dados(&dataset,j);
    }
//Vetor contendo as classes corretas de cada exemplo
    igraph_vector_init(gabarito, l);
    for(i=0 ; i<l ; i++)
            VECTOR(*gabarito)[i] = MATRIX(dataset,i,rotulo_mode);
//Gera a matriz de pesos (l x l)
    igraph_matrix_init(matrix_W, l, l);
    for(i=0 ; i<l ; i++){
        for(j=0 ; j<l ; j++){
            if(i==j)
                MATRIX(*matrix_W, i, j) = 0;
            else
                MATRIX(*matrix_W, i, j) = similarity(&dataset,i,j,ini,fim);
        }
    }
    fclose(instream);
    igraph_matrix_destroy(&dataset);
    return(0);
}

//Aplica a tecnica Knn variando k, a % de pre-rotulados e a quantidade de execucoes
void tecnica_KNN(char title[], int l, int c, int classes, int rotulo_mode, int separador){
    char str[60];
    int k,r,i=0;
    float labeled, taxa_LGC=0, taxa_FH=0, k_med=0;
    srand(time(0));
    igraph_vector_t gabarito;   //rotulos corretos
    igraph_matrix_t matrix_W;   //matriz de similaridade
    igraph_matrix_t matrix_KNN; //matriz de adjacencia apos aplicada a tecnica KNN
    igraph_vector_t rotulos;    //labeled e unlabeled
    igraph_matrix_t dados;      //parametros e dados a serem armazenados
    igraph_matrix_init(&dados, 1, 8);
    gerar_matrizes_fixas(&gabarito,&matrix_W,title,l,c,rotulo_mode,separador);
    for(k=k_inicio ; k<=k_fim ; k++){                                   //grafos com diversos K
        for(labeled=rot_inicio ; labeled<=rot_fim ; labeled+=0.01){     //percentual de dados pre-rotulados
            for(r=0 ; r<exec ; r++){                                    //exec execucoes para cada configuracao acima
                gerar_rede_KNN(&matrix_W,&matrix_KNN,k);
                k_med = grau_medio(&matrix_KNN);
                selecionar_pre_rotulados(&gabarito,&rotulos,labeled,classes);
                taxa_LGC = local_global_consistency(&gabarito,&rotulos,&matrix_KNN,classes,l);
                taxa_FH = Fharmonicas(&gabarito,&rotulos,&matrix_KNN,classes,l);
                armazena_parametros(&dados,(float)k,labeled,i,k_med,taxa_LGC,taxa_FH);
                i++;
            }
        }
    }
    igraph_matrix_remove_row(&dados,i);
    strcpy (str,"Knn-");
    strcat(str,title);
    armazena_dados(str,&dados);
    igraph_matrix_destroy(&matrix_KNN);
    igraph_vector_destroy(&rotulos);
    igraph_matrix_destroy(&matrix_W);
    igraph_vector_destroy(&gabarito);
    igraph_matrix_destroy(&dados);
}

//Aplica o corte epsilon variando e, a % de pre-rotulados e a quantidade de execucoes
void corte_epsilon(char title[], int l, int c, int classes, int rotulo_mode, int separador){
    char str[60];
    int r,i=0;
    float e, labeled, taxa_LGC=0, taxa_FH=0, k_med=0;
    srand(time(0));
    igraph_vector_t gabarito;   //rotulos corretos
    igraph_matrix_t matrix_W;   //matriz de similaridade
    igraph_matrix_t matrix_e;   //matriz apos aplicada o corte epsilon
    igraph_vector_t rotulos;    //labeled e unlabeled
    igraph_matrix_t dados;      //parametros e dados a serem armazenados
    igraph_matrix_init(&dados, 1, 8);
    gerar_matrizes_fixas(&gabarito,&matrix_W,title,l,c,rotulo_mode,separador);
    for(e=e_inicio ; e<=e_fim ; e+=0.01){                               //grafos com diversos e
        for(labeled=rot_inicio ; labeled<=rot_fim ; labeled+=0.01){     //percentual de dados pre-rotulados
            for(r=0 ; r<exec ; r++){                                    //exec execucoes para cada configuracao acima
                gerar_rede_e(&matrix_W,&matrix_e,e);
                k_med = grau_medio(&matrix_e);
                selecionar_pre_rotulados(&gabarito,&rotulos,labeled,classes);
                taxa_LGC = local_global_consistency(&gabarito,&rotulos,&matrix_e,classes,l);
                taxa_FH = Fharmonicas(&gabarito,&rotulos,&matrix_e,classes,l);
                armazena_parametros(&dados,e,labeled,i,k_med,taxa_LGC,taxa_FH);
                i++;
            }
        }
    }
    igraph_matrix_remove_row(&dados,i);
    strcpy (str,"Corte_e-");
    strcat(str,title);
    armazena_dados(str,&dados);
    igraph_matrix_destroy(&matrix_e);
    igraph_vector_destroy(&rotulos);
    igraph_matrix_destroy(&matrix_W);
    igraph_vector_destroy(&gabarito);
    igraph_matrix_destroy(&dados);
}

//                                                                ultima    -> 1 /              virgula -> 1
//nome do dataset / linhas / colunas / qntd. classes / rotulos na 1ª coluna -> 0 / separado por espaço  -> 0
int main(void){
//    tecnica_KNN("Iris_t.txt",9,5,3,1,1);
//    corte_epsilon("Iris_t.txt",9,5,3,1,1);

    tecnica_KNN("Breast_Tissue.txt",106,10,6,0,0);
    corte_epsilon("Breast_Tissue.txt",106,10,6,0,0);

    tecnica_KNN("Iris.txt",150,5,3,1,1);
    corte_epsilon("Iris.txt",150,5,3,1,1);

    tecnica_KNN("Wine.txt",178,14,3,0,0);
    corte_epsilon("Wine.txt",178,14,3,0,0);

    tecnica_KNN("Seeds.txt",210,8,3,1,0);
    corte_epsilon("Seeds.txt",210,8,3,1,0);

    tecnica_KNN("Glass_Identification.txt",214,10,6,1,0);
    corte_epsilon("Glass_Identification.txt",214,10,6,1,0);

    tecnica_KNN("SPECTF_Heart.txt",267,45,2,0,1);
    corte_epsilon("SPECTF_Heart.txt",267,45,2,0,1);

    tecnica_KNN("Haberman_Survival.txt",306,4,2,1,1);
    corte_epsilon("Haberman_Survival.txt",306,4,2,1,1);

    tecnica_KNN("Ecoli.txt",336,8,8,1,0);
    corte_epsilon("Ecoli.txt",336,8,8,1,0);

    tecnica_KNN("Libras_Movement.txt",360,91,15,1,0);
    corte_epsilon("Libras_Movement.txt",360,91,15,1,0);

    tecnica_KNN("Blood_Transfusion_Service_Center.txt",748,5,2,1,0);
    corte_epsilon("Blood_Transfusion_Service_Center.txt",748,5,2,1,0);

    tecnica_KNN("Semeion_Handwritten_Digit.txt",1593,257,10,1,0);
    corte_epsilon("Semeion_Handwritten_Digit.txt",1593,257,10,1,0);

    tecnica_KNN("Wine_Quality_red.txt",1599,12,11,1,0);
    corte_epsilon("Wine_Quality_red.txt",1599,12,11,1,0);

    tecnica_KNN("Image_Segmentation.txt",2310,20,7,0,1);
    corte_epsilon("Image_Segmentation.txt",2310,20,7,0,1);

    tecnica_KNN("Statlog_Image_Segmentation.txt",2310,20,7,1,0);
    corte_epsilon("Statlog_Image_Segmentation.txt",2310,20,7,1,0);

    tecnica_KNN("Wine_Quality_white.txt",4898,12,11,1,0);
    corte_epsilon("Wine_Quality_white.txt",4898,12,11,1,0);

    tecnica_KNN("Optical_Recognition_of_Handwritten_Digits.txt",5620,65,10,1,1);
    corte_epsilon("Optical_Recognition_of_Handwritten_Digits.txt",5620,65,10,1,1);

    tecnica_KNN("Wall-Following_Robot_Navigation_Data.txt",6435,25,4,1,1);
    corte_epsilon("Wall-Following_Robot_Navigation_Data.txt",6435,25,4,1,1);

    tecnica_KNN("Statlog_Landsat_Satellite.txt",6435,37,7,1,0);
    corte_epsilon("Statlog_Landsat_Satellite.txt",6435,37,7,1,0);

    tecnica_KNN("MAGIC_Gamma_Telescope.txt",19020,11,2,1,1);
    corte_epsilon("MAGIC_Gamma_Telescope.txt",19020,11,2,1,1);

    tecnica_KNN("Letter_Recognition.txt",20000,17,26,1,1);
    corte_epsilon("Letter_Recognition.txt",20000,17,26,1,1);

    return(0);
}

