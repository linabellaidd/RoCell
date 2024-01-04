import torch.nn as nn
class AutoEncoder(nn.Module):
    def __init__(self, dims):
        super(AutoEncoder, self).__init__()
        encoder_layers, decoder_layers = [], []
        for i in range(len(dims)):
            if i < len(dims)-1:
                encoder_layers.append(nn.Linear(dims[i], dims[i+1]))
                decoder_layers.append(nn.Linear(dims[len(dims)-1-i], dims[len(dims)-2-i]))
        
        self.encoder_layers = nn.Sequential(*encoder_layers)
        self.decoder_layers = nn.Sequential(*decoder_layers)
        self.leakyrelu = nn.LeakyReLU(0.1)
        self.dropout = nn.Dropout(0.2)
    
    def forward(self, x):
        for encoder_layer in self.encoder_layers:
            x = self.dropout(self.leakyrelu(encoder_layer(x)))
            
        latent_layer = x
        
        for decoder_layer in self.decoder_layers:
            x = self.dropout(self.leakyrelu(decoder_layer(x)))

        return x, latent_layer