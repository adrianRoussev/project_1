Rails.application.routes.draw do
  # Define your application routes per the DSL in https://guides.rubyonrails.org/routing.html
#routes.rb
  
  # Defines the root path route ("/")
  # root "articles#index"

  get "/", to: "welcome#index"

  get "/proteins", to: "proteins#index"
  get "/proteins/new", to: "proteins#new"
  post "/proteins", to: "proteins#create"
  get "/proteins/:id", to: "proteins#show"
  get '/proteins/:id/edit', to: 'proteins#edit'
  patch '/proteins/:id', to: 'proteins#update'
  delete '/proteins/:id', to: 'proteins#destroy'
  get "/inhibitors", to: "inhibitors#index"
  get "/inhibitors/:id", to: "inhibitors#show"
  get "/proteins/:id/inhibitors", to: "protein_inhibitors#index"
  get "/proteins/:id/inhibitors/new", to: "protein_inhibitors#new"
  post "/proteins/:id/inhibitors", to: "protein_inhibitors#create"
  get "/inhibitors/:id/edit", to: "inhibitors#edit"
  patch "/inhibitors/:id", to: "inhibitors#update"
  get "/proteins/:id/inhibitors", to: "protein_inhibitors#index"
  delete "/inhibitors/:id", to: "inhibitors#destroy"
  get "/compounds", to: "compounds#index"
  get "/compounds/new", to: "compounds#new"
  post "/compounds", to: "compounds#create"
  get "/compounds/:id", to: "compounds#show"
  get '/compounds/:id/edit', to: 'compounds#edit'
  patch '/compounds/:id', to: 'compounds#update'
  delete '/compounds/:id', to: 'compounds#destroy'
end