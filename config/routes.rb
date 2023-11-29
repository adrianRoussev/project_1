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
  get "inhibitors/:id/compoundable", to: "inhibitors#new_compoundable"
  post "/inhibitors/:id", to: "inhibitors#create_compoundable"
  get "/inhibitors/:id/reactions", to: "reactions#index"

  get "/proteins/:id/inhibitors", to: "protein_inhibitors#index"
  delete "/inhibitors/:id", to: "inhibitors#destroy"
  get "/compounds", to: "compounds#index"
  get "/compounds/new", to: "compounds#new"
  post "/compounds", to: "compounds#create"
  get "/compounds/:id", to: "compounds#show"
  get '/compounds/:id/edit', to: 'compounds#edit'
  patch '/compounds/:id', to: 'compounds#update'
  delete '/compounds/:id', to: 'compounds#destroy'
  get "/compounds/:id/inhibitor/new", to: "compounds#new_inhibitor"
  post "/compounds/:id/:protein_id/inhibitor", to: "compounds#create_inhibitor"

  get "/reactions", to: "reactions#index"
  get "/reactions/:id", to: "reactions#show"
  get "/inhibitors/:inhibitor_id/reactions/new", to: "reactions#new"
  post "/inhibitors/:inhibitor_id/reactions", to: "reactions#create"
  post "/reactions/predict", to: "reactions#predict"
  get "/reactions/:id", to: "reactions#show"
  get '/reactions/:id/edit', to: 'reactions#edit'
  patch '/reactions/:id', to: 'reactions#update'
  delete '/reactions/:id', to: 'reactions#destroy'

end