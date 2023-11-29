module Compoundable
    extend ActiveSupport::Concern
  
    included do
      has_many :compounds, :as => :compoundable
    end
  end